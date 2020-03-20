
package newance.mzjava.mol.modification.unimod.jaxb;

import javax.xml.bind.annotation.*;


/**
 * <p>Java class for unimod_t complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="unimod_t">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element name="elements" type="{http://www.unimod.org/xmlns/schema/unimod_2}elements_t" minOccurs="0"/>
 *         &lt;element name="modifications" type="{http://www.unimod.org/xmlns/schema/unimod_2}modifications_t" minOccurs="0"/>
 *         &lt;element name="amino_acids" type="{http://www.unimod.org/xmlns/schema/unimod_2}amino_acids_t" minOccurs="0"/>
 *         &lt;element name="mod_bricks" type="{http://www.unimod.org/xmlns/schema/unimod_2}mod_bricks_t" minOccurs="0"/>
 *       &lt;/sequence>
 *       &lt;attribute name="majorVersion" use="required" type="{http://www.w3.org/2001/XMLSchema}unsignedShort" fixed="2" />
 *       &lt;attribute name="minorVersion" use="required" type="{http://www.unimod.org/xmlns/schema/unimod_2}minorVersion_t" />
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "unimod_t", namespace = "http://www.unimod.org/xmlns/schema/unimod_2", propOrder = {
    "modifications"
})
public class UnimodT {

    @XmlElement(namespace = "http://www.unimod.org/xmlns/schema/unimod_2")
    protected ModificationsT modifications;
    @XmlAttribute(name = "majorVersion", required = true)
    @XmlSchemaType(name = "unsignedShort")
    protected int majorVersion;
    @XmlAttribute(name = "minorVersion", required = true)
    protected int minorVersion;

    /**
     * Gets the value of the modifications property.
     * 
     * @return
     *     possible object is
     *     {@link ModificationsT }
     *     
     */
    public ModificationsT getModifications() {
        return modifications;
    }

    /**
     * Sets the value of the modifications property.
     * 
     * @param value
     *     allowed object is
     *     {@link ModificationsT }
     *     
     */
    public void setModifications(ModificationsT value) {
        this.modifications = value;
    }

    /**
     * Gets the value of the majorVersion property.
     * 
     */
    public int getMajorVersion() {
        return majorVersion;
    }

    /**
     * Sets the value of the majorVersion property.
     * 
     */
    public void setMajorVersion(int value) {
        this.majorVersion = value;
    }

    /**
     * Gets the value of the minorVersion property.
     * 
     */
    public int getMinorVersion() {
        return minorVersion;
    }

    /**
     * Sets the value of the minorVersion property.
     * 
     */
    public void setMinorVersion(int value) {
        this.minorVersion = value;
    }

}
