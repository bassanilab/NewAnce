/*
Copyright (C) SIB - Swiss Institute of Bioinformatics, Lausanne, Switzerland
Copyright (C) LICR - Ludwig Institute of Cancer Research, Lausanne, Switzerland
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
*/

package newance.util;

import org.apache.commons.cli.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.InvalidPathException;

/**
* Copyright (C) 2019, SIB/LICR. All rights reserved  *  * SIB, Swiss Institute of Bioinformatics  * Ludwig Institute for Cancer Research (LICR)  *  * Redistribution and use in source and binary forms, with or without  * modification, are permitted provided that the following conditions are met:  * Redistributions of source code must retain the above copyright notice, this  * list of conditions and the following disclaimer. Redistributions in binary  * form must reproduce the above copyright notice, this list of conditions and  * the following disclaimer in the documentation and/or other materials provided  * with the distribution. Neither the name of the SIB/LICR nor the names of  * its contributors may be used to endorse or promote products derived from this  * software without specific prior written permission.  *  * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  * ARE DISCLAIMED. IN NO EVENT SHALL SIB/LICR BE LIABLE FOR ANY DIRECT,  * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES  * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS  * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public abstract class ExecutableOptions {

    protected Options cmdLineOpts;
    protected boolean optionsSet;
    protected String version;

    abstract protected void createOptions();

    protected boolean checkBooleanOption(CommandLine line, String optionId) throws MissingOptionException {

        if (line.hasOption(optionId)) {
            return true;
        }

        return false;
    }

    // if option is present in cmd line, returns option string if option has argument and "true" otherwise.
    // if option is not present in cmd line, returns empty string if option has argument and "false" otherwise.
    protected String getOptionString(CommandLine line, String optionId) throws MissingOptionException{

        String string = "";
        if( line.hasOption( optionId ) ) {

            if (cmdLineOpts.getOption(optionId).hasArg()) {
                string = line.getOptionValue(optionId);

                if (string.isEmpty()) {
                    String optStr = "-" + optionId + ",--" + cmdLineOpts.getOption(optionId).getLongOpt();
                    throw new MissingOptionException(line.getOptionValue(optionId) + " for option " + optStr + ". Missing option string");
                }
            } else {
                string = "true";
            }
        } else {

            if (!cmdLineOpts.getOption(optionId).hasArg()) {
                string = "false";
            }
        }

        return string;
    }

    protected void checkHelpOption(String[] args, String optionId) {

        for (String arg : args) {
            if (arg.equalsIgnoreCase(optionId) || arg.equalsIgnoreCase(cmdLineOpts.getOption(optionId).getLongOpt())) {

                HelpFormatter formatter = new HelpFormatter();
                formatter.setWidth(150);

                for (Option option : cmdLineOpts.getOptions()) {
                    option.setRequired(false);
                }

                System.out.println("");
                System.out.println("*************************************************************************************************************************");
                System.out.println("**                                            NewAnce command line help                                                **");
                System.out.println("*************************************************************************************************************************");
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

    protected boolean checkReadParamsOption(String[] args,  String optionId) {

        for (int i=0;i<args.length-1;i++) {
            if (args[i].equalsIgnoreCase(optionId) || (cmdLineOpts.getOption(optionId)!=null && args[i].equalsIgnoreCase(cmdLineOpts.getOption(optionId).getLongOpt()))) {
                if (!(new File(args[i+1]).exists())) {
                    String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                    throw new InvalidPathException(args[i+1]+" for option "+optStr+". ", "File does not exist");
                }

                NewAnceParams.getInstance().add("readParamsFile",args[i+1]);
                NewAnceParams.getInstance().read(args[i+1]);
                return true;
            }
        }

        return false;
    }

    public ExecutableOptions parseOptions(String[] args) throws ParseException {

        if (!optionsSet) {
            CommandLineParser parser = new DefaultParser();
            CommandLine line = parser.parse(cmdLineOpts, args);
            check(line);
        }

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

    public ExecutableOptions init(String[] args) throws IOException {

        optionsSet = false;
        if (!checkReadParamsOption(args,"-rP")) {
            checkHelpOption(args, "-h");
            checkVersionOption(args, version, "-v");
        } else {
            optionsSet = true;
        }
        return this;
    }

    protected abstract void check(CommandLine line) throws ParseException;
    public abstract int run() throws IOException;

}
