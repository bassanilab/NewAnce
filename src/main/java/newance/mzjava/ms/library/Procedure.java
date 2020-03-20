package newance.mzjava.ms.library;

/**
 * Interface for procedures with one parameter.
 *
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public interface Procedure<T> {

    /**
     * Executes this procedure.
     *
     * @param t the input argument
     */
    void execute(T t);
}
