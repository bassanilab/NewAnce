package newance.mzjava.mol.modification;

import newance.mzjava.mol.Peptide;
import newance.util.NewAnceParams;
import org.junit.Assert;
import org.junit.Test;

/**
 * Created by markusmueller on 10.08.21.
 */
public class NameModificationResolverTest {

    @Test
    public void test_resolve() {
        NameModificationResolver resolver =
                new NameModificationResolver(NewAnceParams.getInstance().getModifications());

        Assert.assertEquals("Phospho", resolver.resolve("Phospho").get().getLabel());
        Assert.assertEquals("Cysteinyl", resolver.resolve("Cysteinyl").get().getLabel());
        Assert.assertEquals("Oxidation", resolver.resolve("Oxidation").get().getLabel());
        Assert.assertEquals("Acetyl", resolver.resolve("Acetyl").get().getLabel());
    }

    @Test
    public void test_resolve2() {
        NameModificationResolver resolver =
                new NameModificationResolver(NewAnceParams.getInstance().getModifications());

        Assert.assertEquals(80, resolver.resolve("Phospho").get().getMolecularMass(), 0.1);
        Assert.assertEquals(119, resolver.resolve("Cysteinyl").get().getMolecularMass(), 0.1);
        Assert.assertEquals(16, resolver.resolve("Oxidation").get().getMolecularMass(), 0.1);
        Assert.assertEquals(42, resolver.resolve("Acetyl").get().getMolecularMass(), 0.1);
        Assert.assertEquals("Cysteinyl", resolver.resolve("Cysteinyl").get().getLabel());
        Assert.assertEquals("Oxidation", resolver.resolve("Oxidation").get().getLabel());
        Assert.assertEquals("Acetyl", resolver.resolve("Acetyl").get().getLabel());
    }

    @Test
    public void test_resolve3() {
        NameModificationResolver resolver =
                new NameModificationResolver(NewAnceParams.getInstance().getModifications());

        Peptide p = Peptide.parse("PEPT(Phospho)IDE", resolver);
        Assert.assertEquals("PEPT(Phospho)IDE", p.toString());
    }

}
