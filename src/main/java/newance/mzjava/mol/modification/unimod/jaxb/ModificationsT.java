
package newance.mzjava.mol.modification.unimod.jaxb;


import newance.mzjava.mol.modification.unimod.UnimodMod;

import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;
import java.util.ArrayList;
import java.util.List;


/**
 * <p>Java class for modifications_t complex type.
 * 
 * <p>The following schema fragment specifies the expected content contained within this class.
 * 
 * <pre>
 * &lt;complexType name="modifications_t">
 *   &lt;complexContent>
 *     &lt;restriction base="{http://www.w3.org/2001/XMLSchema}anyType">
 *       &lt;sequence>
 *         &lt;element name="mod" type="{http://www.unimod.org/xmlns/schema/unimod_2}mod_t" maxOccurs="unbounded" minOccurs="0"/>
 *       &lt;/sequence>
 *     &lt;/restriction>
 *   &lt;/complexContent>
 * &lt;/complexType>
 * </pre>
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "modifications_t", namespace = "http://www.unimod.org/xmlns/schema/unimod_2", propOrder = {
    "mod"
})
public class ModificationsT {

    @XmlElement(namespace = "http://www.unimod.org/xmlns/schema/unimod_2")
    @XmlJavaTypeAdapter(ModificationTypeAdapter.class)
    protected List<UnimodMod> mod;

    /**
     * Gets the value of the mod property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the mod property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getMod().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link ModT }
     * 
     * 
     */
    public List<UnimodMod> getMod() {
        if (mod == null) {
            mod = new ArrayList<UnimodMod>();
        }
        return this.mod;
    }

}
