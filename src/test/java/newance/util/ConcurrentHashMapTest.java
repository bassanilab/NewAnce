package newance.util;

import org.junit.Test;

import java.util.concurrent.ConcurrentHashMap;

/**
 * Copyright (C) 2019
 * @author Markus Müller
 * @Institutions: SIB, Swiss Institute of Bioinformatics; Ludwig Institute for Cancer Research
 */

public class ConcurrentHashMapTest {

    @Test
    public void test_merge() {
        ConcurrentHashMap<Integer, String> conmap = new ConcurrentHashMap<Integer, String>();

        conmap.put(1, "Java");
        conmap.put(2, "php");
        conmap.put(3, ".net");
        conmap.put(5, "python");

        ConcurrentHashMap<Integer, String>  conmap2 = new ConcurrentHashMap<Integer, String>();

        conmap2.put(1, "java");
        conmap2.put(2, "C++");
        conmap2.put(3, "Rubi");
        conmap2.put(6, "Java Script");

        conmap2.forEach(
                (key, value) -> conmap.merge( key, value, (v1, v2) -> v1.equalsIgnoreCase(v2) ? v1 : v1 + "," + v2)
        );

        System.out.println(conmap);

    }
}
