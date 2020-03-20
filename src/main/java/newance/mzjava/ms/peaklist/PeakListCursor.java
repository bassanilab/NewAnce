package newance.mzjava.ms.peaklist;

import java.util.List;

import static com.google.common.base.Preconditions.checkElementIndex;

/**
* @author Oliver Horlacher
* @version sqrt -1
*/
class PeakListCursor<A extends PeakAnnotation> implements PeakCursor<A> {

    private int cursorIndex = -1;

    private final PeakList<A> peakList;
    private final int size;

    PeakListCursor(PeakList<A> peakList) {

        this.peakList = peakList;
        this.size = peakList.size();
    }

    @Override
    public int size() {

        return size;
    }

    @Override
    public boolean isEmpty() {

        return size == 0;
    }

    @Override
    public void resetCursor() {

        cursorIndex = -1;
    }

    @Override
    public boolean next() {

        if (cursorIndex < size - 1) {

            cursorIndex++;

            return true;
        } else {

            return false;
        }
    }

    @Override
    public boolean previous() {

        cursorIndex = Math.max(cursorIndex - 1, -1);
        return cursorIndex > -1;
    }

    @Override
    public List<A> currAnnotations() {

        checkElementIndex(cursorIndex, size);

        return peakList.getAnnotations(cursorIndex);
    }

    @Override
    public boolean canPeek(int n) {

        return cursorIndex + n < size;
    }

    /**
     * Move the cursor to the m/z that is closest to <code>mz</code>
     *
     * @param mz the m/z value to which the cursor is to be moved
     */
    @Override
    public void moveToClosest(double mz) {

        cursorIndex = getClosestIndex(mz);
    }

    @Override
    public void moveBefore(double mz) {

        cursorIndex = getClosestIndex(mz);

        if (currMz() >= mz) {

            previous();
        }
    }

    @Override
    public void movePast(double mz) {

        if (cursorIndex == -1) {

            cursorIndex++;
        }

        while (cursorIndex < size && getMz(cursorIndex) <= mz) {

            cursorIndex++;
        }
    }

    @Override
    public int getClosestIndex(double mz) {

        return peakList.getClosestIndex(mz);
    }

    @Override
    public String toString() {

        return currMz() + ", " + currIntensity();
    }

    @Override
    public boolean next(double mz) {

        if (cursorIndex == -1 || peakList.getMz(cursorIndex) <= mz) {

            cursorIndex++;
        }

        boolean more = cursorIndex < size;

        cursorIndex = Math.min(cursorIndex, size - 1);

        return more;
    }

    @Override
    public double currIntensity() {

        return getIntensity(cursorIndex);
    }

    @Override
    public double currMz() {

        return getMz(cursorIndex);
    }

    @Override
    public double peekMz(int n) {

        return getMz(cursorIndex + n);
    }

    @Override
    public double lastMz() {

        return getMz(size - 1);
    }

    @Override
    public double lastIntensity() {

        return getIntensity(size - 1);
    }

    @Override
    public double peekIntensity(int n) {

        return getIntensity(cursorIndex + n);
    }

    @Override
    public double getMz(int index) {

        checkElementIndex(index, size);

        return peakList.getMz(index);
    }

    @Override
    public double getIntensity(int index) {

        checkElementIndex(index, size);

        return peakList.getIntensity(index);
    }
}
