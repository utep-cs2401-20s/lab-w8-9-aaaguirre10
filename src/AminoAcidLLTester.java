import org.junit.jupiter.api.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.jupiter.api.Assertions.assertArrayEquals;

public class AminoAcidLLTester {


    @Test
    public void aminoAcidCompare(){
        AminoAcidLL a = AminoAcidLL.createFromRNASequence("CAUUUGAUU");
        AminoAcidLL b = AminoAcidLL.createFromRNASequence("CAUUUGAUU");
        a = AminoAcidLL.sort(a);
        b = AminoAcidLL.sort(b);
        assertEquals(0, a.aminoAcidCompare(b));
    }

    @Test
    public void aminoAcidCounts(){
        int[] expected = {2, 2, 1};
        String aminoTester = "GCGGCGCAUCAUUGU";
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence(aminoTester);
        assertArrayEquals(expected, aminoAcid.aminoAcidCounts());
    }

    @Test
    public void aminoAcidCounts2(){
        int[] expected = {1, 1, 1};
        String aminoTester = "GGAGAAGUC";
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence(aminoTester);
        assertArrayEquals(expected, aminoAcid.aminoAcidCounts());
    }

    @Test
    public void aminoAcidList(){
        char[] expected = {'R','E','D','S'};
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("AGGGAGGACUCA");
        assertArrayEquals(expected, aminoAcid.aminoAcidList());
    }

    @Test
    public void aminoAcidList2(){
        char[] expected = {'A','T','C'};
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("GCCUACUGU");
        assertArrayEquals(expected, aminoAcid.aminoAcidList());
    }

    @Test
    public void createFromRNASequence(){
        String expected = "PLA";
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("CCGUUGGCACUGUUG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), aminoAcid.aminoAcid);
            aminoAcid = aminoAcid.next;
        }
    }

    @Test
    public void createFromRNASequence2(){
        String expected = "ATELRS";
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("GCUACGGAGCUUCGGAGCUAG");
        for (int i = 0; i < expected.length(); i++) {
            assertEquals(expected.charAt(i), aminoAcid.aminoAcid);
            aminoAcid = aminoAcid.next;
        }
    }

    @Test
    public void sort(){
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("GAGGAGACC");
        aminoAcid = AminoAcidLL.sort(aminoAcid);
        assertEquals(true, aminoAcid.isSorted());

    }

    @Test
    public void sort2(){
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("AAA");
        aminoAcid = AminoAcidLL.sort(aminoAcid);
        assertEquals(true, aminoAcid.isSorted());

    }

    @Test

    public void isSorted(){
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("GAGGAGACCACCUGCGACUAC");
        assertEquals(false, aminoAcid.isSorted());
    }

    @Test
    public void isSorted2(){
        AminoAcidLL aminoAcid = AminoAcidLL.createFromRNASequence("GCUACGGAGCUUCGGAGCUAG");
        assertEquals(false, aminoAcid.isSorted());
    }

}