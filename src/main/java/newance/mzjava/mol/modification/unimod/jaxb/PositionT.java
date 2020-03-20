
package newance.mzjava.mol.modification.unimod.jaxb;

import javax.xml.bind.annotation.XmlEnum;
import javax.xml.bind.annotation.XmlEnumValue;
import javax.xml.bind.annotation.XmlType;


/**
 * <p>Java class for position_t.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * <p>
 * <pre>
 * &lt;simpleType name="position_t">
 *   &lt;restriction base="{http://www.w3.org/2001/XMLSchema}string">
 *     &lt;enumeration value="Anywhere"/>
 *     &lt;enumeration value="Any N-term"/>
 *     &lt;enumeration value="Any C-term"/>
 *     &lt;enumeration value="Protein N-term"/>
 *     &lt;enumeration value="Protein C-term"/>
 *   &lt;/restriction>
 * &lt;/simpleType>
 * </pre>
 * 
 */
@XmlType(name = "position_t", namespace = "http://www.unimod.org/xmlns/schema/unimod_2")
@XmlEnum
public enum PositionT {

    @XmlEnumValue("Anywhere")
    ANYWHERE("Anywhere"),
    @XmlEnumValue("Any N-term")
    ANY_N_TERM("Any N-term"),
    @XmlEnumValue("Any C-term")
    ANY_C_TERM("Any C-term"),
    @XmlEnumValue("Protein N-term")
    PROTEIN_N_TERM("Protein N-term"),
    @XmlEnumValue("Protein C-term")
    PROTEIN_C_TERM("Protein C-term");
    private final String value;

    PositionT(String v) {
        value = v;
    }

    public String value() {
        return value;
    }

    public static PositionT fromValue(String v) {
        for (PositionT c: PositionT.values()) {
            if (c.value.equals(v)) {
                return c;
            }
        }
        throw new IllegalArgumentException(v);
    }

}
