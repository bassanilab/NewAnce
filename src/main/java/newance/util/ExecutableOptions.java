package newance.util;

import org.apache.commons.cli.*;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.nio.file.InvalidPathException;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

/**
 * Created by markusmueller on 07.12.18.
 */
public abstract class ExecutableOptions {

    protected Options cmdLineOpts;

    abstract protected void createOptions();

    protected boolean checkBooleanOption(CommandLine line, String optionId) throws MissingOptionException {

        if (line.hasOption(optionId)) {
            return true;
        }

        return false;
    }

    protected Set<String> checkStringListOption(CommandLine line, String optionId) throws MissingOptionException {

        Set<String> strSet = new HashSet<>();

        String string = "";
        if (line.hasOption(optionId)) {
            string = line.getOptionValue(optionId);

            if (string.isEmpty()) {
                String optStr = "-" + optionId + ",--" + cmdLineOpts.getOption(optionId).getLongOpt();
                throw new MissingOptionException(line.getOptionValue(optionId) + " for option " + optStr + ". Missing option string");
            }

            String[] items = string.split(",");

            for (String str : items) {
                strSet.add(str);
            }
        }

        return strSet;
    }


    protected String checkStringOption(CommandLine line, String optionId) throws MissingOptionException{

        String string = "";
        if( line.hasOption( optionId ) ) {
            string = line.getOptionValue( optionId );

            if (string.isEmpty()) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new MissingOptionException(line.getOptionValue(optionId)+" for option "+optStr+". Missing option string");
            }
        }

        return string;
    }


    protected String checkDefinedStringOption(CommandLine line, String optionId, Set<String> allowedStrings) throws MissingOptionException{

        String string = "";
        if( line.hasOption( optionId ) ) {
            string = line.getOptionValue( optionId ).toLowerCase();

            if (string.isEmpty()) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new MissingOptionException(line.getOptionValue(optionId)+" for option "+optStr+". Missing option string");
            } else if (!allowedStrings.contains(string)) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new MissingOptionException(line.getOptionValue(optionId)+" for option "+optStr+". Option string has to be one of "+allowedStrings.toString());
            }
        }

        return string;
    }

    protected int checkIntOption(CommandLine line,String optionId, int min, int max, int defaultValue) {

        int value = defaultValue;

        if( line.hasOption( optionId ) ) {

            String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
            try {
                value = Integer.parseInt(line.getOptionValue(optionId));
            } catch(NumberFormatException e) {
                throw new NumberFormatException(optStr+": bad format for integer option value "+line.getOptionValue(optionId));
            }

            if (value<min || value>max) {
                throw new NumberFormatException(optStr+" option value must be >="+min+" && <="+max);
            }
        }

        return value;
    }

    protected float checkFloatOption(CommandLine line,String optionId, float min, float max, float defaultValue) {
        float value = defaultValue;

        if( line.hasOption( optionId ) ) {

            String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
            try {
                value = Float.parseFloat(line.getOptionValue(optionId));
            } catch(NumberFormatException e) {
                throw new NumberFormatException(optStr+": bad format for integer option value "+line.getOptionValue(optionId));
            }

            if (value<min || value>max) {
                throw new NumberFormatException(optStr+" option value must be >="+min+" && <="+max);
            }
        }

        return value;
    }

    protected double checkDoubleOption(CommandLine line,String optionId, double min, double max, double defaultValue) {
        double value = defaultValue;

        if( line.hasOption( optionId ) ) {

            String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
            try {
                value = Double.parseDouble(line.getOptionValue(optionId));
            } catch(NumberFormatException e) {
                throw new NumberFormatException(optStr+": bad format for integer option value "+line.getOptionValue(optionId));
            }

            if (value<min || value>max) {
                throw new NumberFormatException(optStr+" option value must be >="+min+" && <="+max);
            }
        }

        return value;
    }

    protected String checkExistFileOption(CommandLine line, String optionId) {

        String fileName = "";
        if( line.hasOption( optionId ) ) {
            fileName = line.getOptionValue( optionId );

            if (!(new File(fileName).exists())) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+". ", "File does not exist");
            }
        }

        return fileName;
    }

    protected String checkCanWriteFileOption(CommandLine line, String optionId) {

        String fileName = "";
        if( line.hasOption( optionId ) ) {
            fileName = line.getOptionValue( optionId );

            File file = new File(fileName);
            if (file.exists()) {
                if (!file.canWrite()) {
                    String optStr = "-" + optionId + ",--" + cmdLineOpts.getOption(optionId).getLongOpt();
                    throw new InvalidPathException(line.getOptionValue(optionId) + " for option " + optStr+". ", "File does not exist");
                }
            } else {
                if (!(new File(file.getParent())).canWrite()) {
                    String optStr = "-" + optionId + ",--" + cmdLineOpts.getOption(optionId).getLongOpt();
                    throw new InvalidPathException(line.getOptionValue(optionId) + " for option " + optStr+". ", "File does not exist");
                }
            }
        }

        return fileName;
    }

    protected String checkExistDirOption(CommandLine line, String optionId) {

        String fileName = "";
        if( line.hasOption( optionId ) ) {
            fileName = line.getOptionValue( optionId );

            File dir = new File(fileName);
            if (!dir.isDirectory()) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+"."," Directory does not exist. Please create directory.");
            }

        }

        return fileName;
    }

    protected String checkNotExistDirOption(CommandLine line, String optionId) {

        String fileName = "";
        if( line.hasOption( optionId ) ) {
            fileName = line.getOptionValue( optionId );

            File dir = new File(fileName);
            if (dir.isDirectory()) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+"."," Directory does already exist. Please choose another directory or delete old one.");
            }

        }

        return fileName;
    }

    protected String checkDeleteOldDirOption(CommandLine line, String optionId, boolean forceDeleteDir) {

        String dirName = "";
        if( line.hasOption( optionId ) ) {
            dirName = line.getOptionValue( optionId );

            File dir = new File(dirName);
            try {
                if (dir.isDirectory()) {
                    if (forceDeleteDir) {
                        FileUtils.deleteDirectory(dir);
                    } else {
                        String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                        throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+"."," Directory does already exist. Please choose another directory, delete old one or use -fd,--forceDelete option.");
                    }
                }
            } catch (IOException e) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+"."," Cannot writeHistograms hadoop file to directory");
            }
        }

        return dirName;
    }

    protected String checkOutputDirOption(CommandLine line, String optionId) {

        String dirName = "";
        if( line.hasOption( optionId ) ) {
            dirName = line.getOptionValue( optionId );

            File dir = new File(dirName);
            try {
                if (!dir.exists())
                    dir.mkdirs();
            } catch (SecurityException e) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new InvalidPathException(line.getOptionValue(optionId)+" for option "+optStr+"."," Cannot delete directory");
            }
        }

        return dirName;
    }


    protected void checkHelpOption(String[] args, String optionId) {

        for (String arg : args) {
            if (arg.equalsIgnoreCase(optionId) || arg.equalsIgnoreCase(cmdLineOpts.getOption(optionId).getLongOpt())) {

                HelpFormatter formatter = new HelpFormatter();

                for (Option option : cmdLineOpts.getOptions()) {
                    option.setRequired(false);
                }

                System.out.println("*******************************");
                System.out.println("** NewAnce command line help **");
                System.out.println("*******************************");
                System.out.println("");
                System.out.println("");
                formatter.printHelp( this.getClass().getName(), cmdLineOpts );

                System.exit(1);
            }
        }

    }

    protected void checkVersionOption(String[] args, String version, String optionId) {

        for (String arg : args) {
            if (arg.equalsIgnoreCase(optionId) || arg.equalsIgnoreCase(cmdLineOpts.getOption(optionId).getLongOpt())) {

                System.out.println("NewAnce "+version+" \n");

                System.exit(1);
            }
        }

    }

    protected Pattern checkRegExOption(CommandLine line, String optionId) {

        Pattern regExp = null;
        if( line.hasOption( optionId ) ) {
            String regexStr = line.getOptionValue( optionId );

            try {
                if (!regexStr.isEmpty()) regExp = Pattern.compile(regexStr);
            } catch (PatternSyntaxException exception) {
                String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                throw new PatternSyntaxException(regexStr+" invalid regex for option "+optStr,regexStr,exception.getIndex());
            }
        }

        return regExp;
    }


    public ExecutableOptions parseOptions(String[] args) throws ParseException {

        CommandLineParser parser = new DefaultParser();
        CommandLine line = parser.parse(cmdLineOpts, args);

        check(line);

        return this;
    }

    public void printOptions(String[] args, String msg) {
        String cmdl = "";
        for (int i=0;i<args.length;i++)  cmdl += args[i]+" ";
        System.out.println(cmdl);
        System.out.println(msg);
        System.out.println();
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp( this.getClass().getName(), cmdLineOpts );
    }

    protected abstract void check(CommandLine line) throws ParseException;
    public abstract ExecutableOptions init() throws IOException;
    public abstract int run() throws IOException;

}
