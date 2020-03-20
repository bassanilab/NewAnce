/**
 * Copyright (c) 2010, SIB. All rights reserved.
 *
 * SIB (Swiss Institute of Bioinformatics) - http://www.isb-sib.ch Host -
 * http://mzjava.expasy.org
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer. Redistributions in binary
 * form must reproduce the above copyright notice, this list of conditions and
 * the following disclaimer in the documentation and/or other materials provided
 * with the distribution. Neither the name of the SIB/GENEBIO nor the names of
 * its contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL SIB/GENEBIO BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package newance.mzjava.ms.io.mgf;

import gnu.trove.list.array.TIntArrayList;
import newance.mzjava.ms.peaklist.Peak;
import newance.mzjava.ms.peaklist.PeakAnnotation;
import newance.mzjava.ms.peaklist.PeakList;
import newance.mzjava.ms.peaklist.PeakProcessorChain;
import newance.mzjava.ms.spectrum.*;

import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Reader;
import java.net.URI;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * To extend how the TITLE line is handled add more TitleParser instances using addTitleParser
 *
 * @author Oliver Horlacher
 * @author fnikitin
 * @version 1.0
 */
public abstract class AbstractMgfReader<A extends PeakAnnotation, S extends Spectrum<A>> extends AbstractMsReader<A, S> {

    private static final int READ_AHEAD_LIMIT = 4096; //Half of the default buffer size
    private static final String BEGIN_IONS = "BEGIN IONS";

    private static final Pattern PAT_CHARGE = Pattern.compile("([-+]?)(\\d+)([-+]?)");
    private static final Pattern PAT_COMMENT = Pattern.compile("^[#/!;].*");
    private static final Pattern PAT_PEP_MASS = Pattern.compile("^(\\d+\\.?\\d*)\\s*(\\d+\\.?\\d*)?\\s*([\\d\\s*,?\\s*]*)?");
    private static final Pattern PAT_PEP_MASS_CHARGE = Pattern.compile("\\d+");
    private static final Pattern PAT_SCANS = Pattern.compile("\\s*(\\d+-?\\d*)[,\\s]*");
    private static final Pattern PAT_RT_TAG = Pattern.compile("(\\d+\\.?\\d*-?\\d*\\.?\\d*)");
    private static final Pattern PAT_PEAK =
            Pattern.compile("(" + RegexConstants.REAL + "|"
                    + RegexConstants.INTEGER + ")\\s+(" + RegexConstants.REAL + "|"
                    + RegexConstants.INTEGER + ").*");

    //Set to a non null value by parseHeader
    private int[] defaultCharge;

    public AbstractMgfReader(Reader reader, final URI spectraSource, PeakList.Precision precision, PeakProcessorChain<A> processorChain) throws IOException {

        super(reader, spectraSource, precision, processorChain);
    }

    public String getExtension() {

        return "mgf";
    }

    @Override
    protected String parseHeader(ParseContext context) throws IOException {

        LineNumberReader reader = context.getCurrentReader();
        StringBuilder buff = new StringBuilder();
        final String lineSeparator = System.getProperty("line.separator");

        boolean defaultChargeSet = false;

        reader.mark(READ_AHEAD_LIMIT);
        String line = reader.readLine();
        while (line != null && (line.length() == 0 || !Character.isDigit(line.charAt(0))) && !line.startsWith(BEGIN_IONS)) {

            buff.append(line);
            buff.append(lineSeparator);

            if (line.startsWith("CHARGE=")) {

                TIntArrayList charges = extractCharge(line.substring(7));

                defaultCharge = charges.toArray();
                defaultChargeSet = true;
            }

            reader.mark(READ_AHEAD_LIMIT);
            line = reader.readLine();
        }

        if (!defaultChargeSet) defaultCharge = new int[0];

        reader.reset();
        return buff.toString();
    }

    /**
     * Parse all the lines that belong to one MS1 spectrum and the metadata and peaks to the spectrum.
     *
     * @param spectrum the spectrum
     * @param context the parse context
     * @throws IOException if there are IO exceptions
     */
    protected void parseNextMs1Entry(S spectrum, ParseContext context) throws IOException {

        LineNumberReader reader = context.getCurrentReader();

        boolean endTag = false;

        for (String line = reader.readLine(); line != null; line = reader.readLine()) {

            reader.mark(READ_AHEAD_LIMIT);

            if (line.length() > 0 && Character.isDigit(line.charAt(0))) {

                handlePeakLine(spectrum, line, context);
            } else if (line.startsWith(BEGIN_IONS)) {

                // reset to the previous mark (for next call of parseNextEntry() that need a valid entry)
                reader.reset();

                endTag = true;
                break;
            } else if (line.length() == 0) {

                handleEndIonLine(spectrum, line);

                endTag = true;
                break;
            } else if (matches(line, PAT_COMMENT)) {

                handleComment(spectrum, line);
            } else {

                handleUnknownLine(spectrum, line);
            }
        }

        if (!endTag) throw new IOException("mandatory 'END IONS' tag is missing!");
    }

    /**
     * Parse all lines that correspond to the next MS2 entry and the metadata and peaks to the spectrum.
     *
     * @param spectrum the spectrum
     * @param context the parse context
     * @throws IOException if there are IO exceptions
     */
    protected void parseNextMs2Entry(S spectrum, ParseContext context) throws IOException {

        LineNumberReader reader = context.getCurrentReader();

        boolean endTag = false;

        String line;

        for (line = reader.readLine(); line != null; line = reader.readLine()) {

            if (line.length() > 0 && Character.isDigit(line.charAt(0))) {

                handlePeakLine(spectrum, line, context);
            } else if (line.startsWith(BEGIN_IONS)) {

                // begin tag has already been met! -> missing END IONS
                // reset to the previous mark (for next call of parseNextEntry() that need a valid entry)
                reader.reset();
                break;
            } else if (line.startsWith("END IONS")) {

                handleEndIonLine(spectrum, line);

                endTag = true;
                break;
            } else if (line.contains("=")) {

                final int index = line.indexOf('=');
                handleTagLine(line.substring(0, index), line.substring(index + 1, line.length()), spectrum);
            } else if (matches(line, PAT_COMMENT)) {

                handleComment(spectrum, line);
            } else {

                handleUnknownLine(spectrum, line);
            }

            reader.mark(READ_AHEAD_LIMIT);
        }

        if (!endTag) throw new IOException("mandatory 'END IONS' tag is missing!");
    }

    /**
     * Allows subclasses to handle lines that occur between BEGIN IONS and END IONS tags.
     * Default implementation does nothing
     *
     * @param line the line
     */
    protected void handleIntermediaryLine(@SuppressWarnings("UnusedParameters") String line) {

    }

    /**
     * Handle intermediary lines between BEGIN IONS and END IONS.
     *
     * @param line the intermediary line to handle
     * @return true if this line must be skipped
     */
    protected boolean skipIntermediaryLine(String line) {

        if(line.isEmpty()) {

            return true;
        } else if (Character.isDigit(line.charAt(0))) { // next MS1 entry: stop skipping!
            return false;
        } else if (line.startsWith(BEGIN_IONS)) { // next MSn entry: stop skipping!
            return false;
        } else {
            handleIntermediaryLine(line);
            return true;
        }
    }

    @Override
    protected S parseNextEntry(ParseContext context) throws IOException {

        LineNumberReader reader = context.getCurrentReader();

        String line;

        // skipping intermediary lines
        for (line = reader.readLine(); line != null; line = reader.readLine()) {

            if (!skipIntermediaryLine(line)) break;

            reader.mark(READ_AHEAD_LIMIT);
        }

        // stop if EOF
        if (line == null) return null;

        // prepare next spectrum
        S spectrum = newSpectrum(context, precision);
        spectrum.getPrecursor().setMzAndCharge(0, defaultCharge);
        spectrum.setMsLevel(1);

        if (line.startsWith(BEGIN_IONS)) {

            handleBeginIonLine(spectrum, line);
            parseNextMs2Entry(spectrum, context);
        } else if (context.getNumberOfParsedEntry() == 0) {

            reader.reset();
            parseNextMs1Entry(spectrum, context);
        }

        return spectrum;
    }

    /**
     * Create a new spectrum. This method allows subclasses to specify what type of spectra is read.
     *
     * @param context the parse context
     * @param precision the precision for the new spectrum
     * @return the newly created spectrum
     */
    protected abstract S newSpectrum(ParseContext context, PeakList.Precision precision);

    /**
     * Convenience method to check if a string matches a pattern
     *
     * @param string  the string
     * @param pattern the patter
     * @return true if the string matches the patter, false otherwise
     */
    private static boolean matches(final String string, final Pattern pattern) {

        Matcher matcher = pattern.matcher(string);
        return matcher.matches();
    }

    /**
     * Handle comment line.  Comment lies start with #, ;, ! or /.
     * Default implementation does nothing.
     *
     * @param spectrum the spectrum currently being read
     * @param line     the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handleComment(@SuppressWarnings("UnusedParameters") S spectrum, @SuppressWarnings("UnusedParameters") String line) {

        return false;
    }

    /**
     * Handle Begin Ion line
     *
     * @param spectrum the spectrum currently being read
     * @param line     the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handleBeginIonLine(S spectrum, @SuppressWarnings("UnusedParameters") String line) {

        spectrum.setMsLevel(2);
        return true;
    }

    /**
     * Handle tag value line.  The tag and value is separated by =
     *
     * @param spectrum the spectrum currently being read
     * @param tag      the string from the start of the lien to the = character
     * @param value    the string from the character after the = to the end of the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handleTagLine(String tag, String value, S spectrum) {

        if (tag.startsWith("TITLE")) {

            return parseTitleTag(value, spectrum);
        } else if (tag.startsWith("PEPMASS")) {

            return parsePepMassTag(value, spectrum);
        } else if (tag.startsWith("CHARGE")) {

            return parseChargeTag(value, spectrum);
        } else if (tag.startsWith("SCANS")) {

            return parseScanTag(value, spectrum);
        } else if (tag.startsWith("RTINSECONDS")) {

            return parseRetentionTimeTag(value, spectrum);
        } else {

            return parseUnknownTag(tag, value, spectrum);
        }
    }

    /**
     * Handle peak line, which is any line that starts with a digit
     *
     * @param spectrum the spectrum currently being read
     * @param line     the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handlePeakLine(S spectrum, String line, ParseContext context) {

        Matcher peakMatcher = PAT_PEAK.matcher(line);
        if (peakMatcher.matches()) {

            double mz = Double.parseDouble(peakMatcher.group(1));
            double intensity = Double.parseDouble(peakMatcher.group(2));

            addPeakToSpectrum(spectrum, mz, intensity, context);
            return true;
        } else {

            return false;
        }
    }

    /**
     * Handle end ion line, which is a line that starts with "END IONS".
     * Default implementation does nothing
     *
     * @param spectrum the spectrum currently being read
     * @param line     the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handleEndIonLine(@SuppressWarnings("UnusedParameters") S spectrum, @SuppressWarnings("UnusedParameters") String line) {

        return false;
    }

    /**
     * Handle line that is not a comment, begin ion, tag, peak or end ion.
     * Default implementation does nothing
     *
     * @param spectrum the spectrum currently being read
     * @param line     the line
     * @return true if the information in the line was processed, false otherwise
     */
    protected boolean handleUnknownLine(@SuppressWarnings("UnusedParameters") S spectrum, @SuppressWarnings("UnusedParameters") String line) {

        return false;
    }

    /**
     * Extracts charges from a string. Any substring that matches [+-]\d+[+-] is said to be a charge
     * the sign of the charge is determined by the sign before or after the digit.  If there is a sign
     * before and after the digit the sign after the digit takes precedent.
     *
     * @param text the String from which to extract the charges
     * @return TIntArrayList containing the charges
     */
    TIntArrayList extractCharge(String text) {

        TIntArrayList charges = new TIntArrayList();

        Pattern pattern = PAT_CHARGE;

        Matcher matcher = pattern.matcher(text);
        while (matcher.find()) {

            int number = Integer.parseInt(matcher.group(2));
            String sign = matcher.group(3);
            if (sign == null || sign.length() == 0) sign = matcher.group(1);

            if (sign != null && sign.length() == 1 && sign.charAt(0) == '-')
                number = number * -1;

            charges.add(number);
        }
        return charges;
    }

    /**
     * Parse the line that starts with a TITLE tag.  This uses the list of title parsers
     * to find a parser that can parse the title. Note that the title parsed value is assign to
     * the spectrum metadata comment by default.
     *
     * @param value    the title value
     * @param spectrum the spectrum being parsed
     * @return true if any of the registered title parses extracted information from the title
     */
    protected abstract boolean parseTitleTag(String value, S spectrum);

    /**
     * Called if there are any tags that are not handled by MgfReader
     * <p/>
     * This is meant to be over ridden by subclasses that want to
     * parse tags that the MgfReader does not know about
     *
     * @param tag      the characters in the line before the first =
     * @param value    the characters in the line after the first =
     * @param spectrum the spectrum that is being read
     * @return true of this handled the tag, false otherwise
     */
    @SuppressWarnings("UnusedParameters")
    //This is meant to be over ridden by subclasses that want to parse tags that the MgfReader does not know about
    protected boolean parseUnknownTag(String tag, String value, S spectrum) {

        return false;
    }

    /**
     * Parse the line that starts with RTINSECONDS
     *
     * @param value    the value
     * @param spectrum the spectrum being parsed
     * @return true if the tag was handled, false otherwise
     */
    protected boolean parseRetentionTimeTag(String value, S spectrum) {

        RetentionTimeList retentionTimeList = new RetentionTimeList();
        Matcher matcher = PAT_RT_TAG.matcher(value);

        while (matcher.find()) {

            String rt = matcher.group(1);

            final int index = rt.indexOf('-');
            if (index != -1) {

                retentionTimeList.add(new RetentionTimeInterval(
                        Double.parseDouble(rt.substring(0, index)),
                        Double.parseDouble(rt.substring(index + 1)),
                        TimeUnit.SECOND
                ));
            } else {

                retentionTimeList.add(new RetentionTimeDiscrete(Double.parseDouble(rt), TimeUnit.SECOND));
            }
        }

        if (!retentionTimeList.isEmpty()) {

            setRetentionTimes(spectrum, retentionTimeList);
            return true;
        } else {

            return false;
        }
    }

    /**
     * Add the retention times contained in the retentionTimeList to the spectrum.
     *
     * @param spectrum the spectrum that is currently being read
     * @param retentionTimeList the retention times that were read
     */
    protected abstract void setRetentionTimes(S spectrum, RetentionTimeList retentionTimeList);

    /**
     * Parse the line that starts with CHARGE
     *
     * @param value    the value
     * @param spectrum the spectrum being parsed
     * @return true if the tag was handled, false otherwise
     */
    protected boolean parseScanTag(String value, S spectrum) {

        ScanNumberList scanNumbers = new ScanNumberList();
        Matcher matcher = PAT_SCANS.matcher(value);

        while (matcher.find()) {

            String scanNumber = matcher.group(1);

            final int index = scanNumber.indexOf('-');
            if (index != -1) {

                scanNumbers.add(
                        Integer.parseInt(scanNumber.substring(0, index)), Integer.parseInt(scanNumber.substring(index + 1))
                );
            } else {

                scanNumbers.add(Integer.parseInt(scanNumber));
            }
        }

        if (!scanNumbers.isEmpty()) {

            setScanNumbers(spectrum, scanNumbers);
            return true;
        } else {

            return false;
        }
    }

    /**
     * Add the scan numbers contained in the scanNumbers to the spectrum.
     *
     * @param spectrum the spectrum that is currently being read
     * @param scanNumbers the scan numbers that were read
     */
    protected abstract void setScanNumbers(S spectrum, ScanNumberList scanNumbers);

    /**
     * Parse the line that starts with CHARGE
     *
     * @param value    the value
     * @param spectrum the spectrum being parsed
     * @return true if the tag was handled, false otherwise
     */
    protected boolean parseChargeTag(String value, S spectrum) {

        final Peak precursor = spectrum.getPrecursor();

        TIntArrayList charges = extractCharge(value);

        precursor.setMzAndCharge(precursor.getMz(), charges.toArray());

        return true;
    }

    /**
     * Parse the line that starts with PEPMASS
     *
     * @param value    the value
     * @param spectrum the spectrum being parsed
     * @return true if the tag was handled, false otherwise
     */
    protected boolean parsePepMassTag(String value, S spectrum) {

        Matcher matcher = PAT_PEP_MASS.matcher(value);
        if (matcher.find()) {

            String mzGroup = matcher.group(1);
            String intensityGroup = matcher.group(2);
            String chargeGroup = matcher.group(3);

            Peak precursor = spectrum.getPrecursor();
            double mz = Double.parseDouble(mzGroup);
            int[] charge = precursor.getChargeList();
            if (intensityGroup != null) precursor.setIntensity(Double.parseDouble(intensityGroup));
            if (chargeGroup != null && chargeGroup.length() > 0) {

                TIntArrayList charges = new TIntArrayList();

                matcher = PAT_PEP_MASS_CHARGE.matcher(chargeGroup);
                while (matcher.find()) {

                    charges.add(Integer.parseInt(matcher.group(0)));
                }

                charge = charges.toArray();
            }
            precursor.setMzAndCharge(mz, charge);
            return true;
        } else {

            return false;
        }
    }
}
