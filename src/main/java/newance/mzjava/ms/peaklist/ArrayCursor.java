package newance.mzjava.ms.peaklist;

import java.util.*;

import static com.google.common.base.Preconditions.*;

/**
 * A cursor for traversing an mz and intensity array
 *
 * @author fnikitin
 * Date: 3/15/13
 */
public class ArrayCursor<A extends PeakAnnotation> implements PeakCursor<A> {

    protected final double[] mzList;
    protected final double[] intensityList;
    protected final int size;

    protected Map<Integer, List<A>> annotationMap;

    protected int cursorIndex = -1;

    public ArrayCursor(double[] mzList, double[] intensityList, int size) {

        this(mzList, intensityList, new HashMap<Integer, List<A>>(), size);
    }

    public ArrayCursor(double[] mzList, double[] intensityList, Map<Integer, List<A>> annotationMap, int size) {

        checkNotNull(mzList);
        checkNotNull(intensityList);
        checkNotNull(annotationMap);
        checkArgument(size <= mzList.length && size <= intensityList.length,
                "size is smaller than arrays");

        this.mzList = Arrays.copyOf(mzList, size);
        this.intensityList = Arrays.copyOf(intensityList, size);
        this.size = size;
        this.annotationMap = annotationMap;
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

        if (annotationMap.containsKey(cursorIndex)) {

            return Collections.unmodifiableList(annotationMap.get(cursorIndex));
        } else {

            return Collections.emptyList();
        }
    }

    @Override
    public boolean canPeek(int n) {

        return cursorIndex + n < size;
    }

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

        cursorIndex = getClosestIndex(mz);

        if (currMz() <= mz) {

            cursorIndex++;
            cursorIndex = Math.min(cursorIndex, size - 1);
        }
    }

    @Override
    public int getClosestIndex(double mz) {

        if (size == 0) {

            return -1;
        }

        int index = Arrays.binarySearch(mzList, 0, size, mz);

        if (index < 0) index = -1 * (index + 1);

        if (index == size) {

            return index - 1;
        } else if (index == 0) {

            return index;
        }

        double ds = mz - getMz(index - 1);
        double dl = getMz(index) - mz;

        return ds < dl ? index - 1 : index;
    }

    @Override
    public boolean next(double mz) {

        if (cursorIndex == -1 || mzList[cursorIndex] <= mz) {

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

        return mzList[index];
    }

    @Override
    public double getIntensity(int index) {

        checkElementIndex(index, size);

        return intensityList[index];
    }

    @Override
    public String toString() {

        return "{"+currMz() + ", " + currIntensity()+", "+currAnnotations()+"}";
    }
}