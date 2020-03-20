package newance.mzjava.mol.modification;

import java.util.AbstractList;

/**
 * @author Oliver Horlacher
 * @version sqrt -1
 */
public class ModificationLists {

    public static final ModificationList EMPTY_MOD_LIST = new EmptyModList();

    private ModificationLists() {}

    private static class EmptyModList extends AbstractList<Modification> implements ModificationList {

        @Override
        public Modification get(int index) {

            throw new IndexOutOfBoundsException("Index: " + index + ", Size: 0");
        }

        @Override
        public int size() {

            return 0;
        }

        @Override
        public double getMolecularMass() {

            return 0;
        }
    }
}
