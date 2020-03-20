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


/**
 * A list of static regular expressions.
 *
 * @author nikitin
 * @version 1.0
 */
public final class RegexConstants {

    /**
     * algebraic sign
     */
    private static final String SIGN = "[+-]";

    /**
     * integer expression
     */
    public static final String INTEGER = SIGN + "?\\d+";

    /**
     * float expression
     */
    public static final String REAL = SIGN + "?\\d+\\.?\\d*(?:[eE][-+]?\\d+)?";

    /**
     * About new lines, sources at http://en.wikipedia.org/wiki/Newline:
     * <p/>
     * Systems based on ASCII or a compatible character set use either LF (Line
     * feed, '\n', 0x0A, 10 in decimal) or CR (Carriage return, '\r', 0x0D, 13
     * in decimal) individually, or CR followed by LF (CR+LF, 0x0D 0x0A). These
     * characters are based on printer commands: The line feed indicated that
     * one line of paper should feed out of the printer, and a carriage return
     * indicated that the printer carriage should return to the beginning of the
     * current line.
     * <p/>
     * <ul>
     * <li>LF: Multics, Unix and Unix-like systems (GNU/Linux, AIX, Xenix, Mac
     * OS X, FreeBSD, etc.), BeOS, Amiga, RISC OS, and others</li>
     * <li>CR+LF: DEC RT-11 and most other early non-Unix, non-IBM OSes, CP/M,
     * MP/M, MS-DOS, OS/2, Microsoft Windows, Symbian OS</li>
     * <li>CR: Commodore 8-bit machines, TRS-80, Apple II family, Mac OS up to
     * version 9 and OS-9</li>
     * </ul>
     */
    public static final String LINE_DELIMITOR = "(?:\r\n|\n|\r)";

    private RegexConstants() {

    }
}
