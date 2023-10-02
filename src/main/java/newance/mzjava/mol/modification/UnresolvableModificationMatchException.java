package newance.mzjava.mol.modification;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class UnresolvableModificationMatchException extends RuntimeException {

    private final ModificationMatch match;

    public UnresolvableModificationMatchException(ModificationMatch match) {

        super("ModificationMatch " + match.toString() + " does not map a known modification");

        this.match = match;
    }

    public ModificationMatch getMatch() {

        return match;
    }
}
