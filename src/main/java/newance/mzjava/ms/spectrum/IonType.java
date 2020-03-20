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
package newance.mzjava.ms.spectrum;

/**
 *
 * @author Oliver Horlacher
 */
public enum IonType {

    a(FragmentType.FORWARD), b(FragmentType.FORWARD), c(FragmentType.FORWARD), d(FragmentType.FORWARD),
    w(FragmentType.REVERSE), x(FragmentType.REVERSE), y(FragmentType.REVERSE), z(FragmentType.REVERSE),
    //Internal ions
    aw_internal(FragmentType.INTERNAL), ax_internal(FragmentType.INTERNAL), ay_internal(FragmentType.INTERNAL), az_internal(FragmentType.INTERNAL),
    bw_internal(FragmentType.INTERNAL), bx_internal(FragmentType.INTERNAL), by_internal(FragmentType.INTERNAL), bz_internal(FragmentType.INTERNAL),
    cw_internal(FragmentType.INTERNAL), cx_internal(FragmentType.INTERNAL), cy_internal(FragmentType.INTERNAL), cz_internal(FragmentType.INTERNAL),
    dw_internal(FragmentType.INTERNAL), dx_internal(FragmentType.INTERNAL), dy_internal(FragmentType.INTERNAL), dz_internal(FragmentType.INTERNAL),
    i(FragmentType.MONOMER), //Immonium
    p(FragmentType.INTACT), //Precursor
    unknown(FragmentType.UNKNOWN);

    private FragmentType fragmentType;


    IonType(FragmentType fragmentType) {

        this.fragmentType = fragmentType;
    }

    public FragmentType getFragmentType() {

        return fragmentType;
    }

    public IonType getComplement() {

        switch (this) {

            case a:
                return x;
            case b:
                return y;
            case c:
                return z;
            case x:
                return a;
            case y:
                return b;
            case z:
                return c;
            default:
                throw new IllegalStateException("Cannot get complement of " + this);
        }
    }
}

