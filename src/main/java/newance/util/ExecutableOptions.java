package newance.util;

import newance.psmcombiner.CometMaxQuantScoreCombiner;
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
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public abstract class ExecutableOptions {

    protected Options cmdLineOpts;
    protected boolean optionsSet;

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

    protected boolean checkReadParamsOption(String[] args,  String optionId) {

        for (int i=0;i<args.length-1;i++) {
            if (args[i].equalsIgnoreCase(optionId) || args[i].equalsIgnoreCase(cmdLineOpts.getOption(optionId).getLongOpt())) {
                if (!(new File(args[i+1]).exists())) {
                    String optStr = "-"+optionId+",--"+cmdLineOpts.getOption(optionId).getLongOpt();
                    throw new InvalidPathException(args[i+1]+" for option "+optStr+". ", "File does not exist");
                }

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
            checkVersionOption(args, NewAnceParams.getInstance().getVersion(), "-v");
        } else {
            optionsSet = true;
        }
        return this;
    }

    protected abstract void check(CommandLine line) throws ParseException;
    public abstract int run() throws IOException;

}
