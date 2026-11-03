/*
 * MORTAR - MOlecule fRagmenTAtion fRamework
 * Copyright (C) 2026  Felix Baensch, Jonas Schaub (felix.j.baensch@gmail.com, jonas.schaub@uni-jena.de)
 *
 * Source code is available at <https://github.com/FelixBaensch/MORTAR>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.mortar.model.fragmentation.algorithm;

import de.unijena.cheminf.mortar.model.util.ChemUtil;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Generator of a systematic set of hydrocarbon test molecules for the {@link AlkylStructureFragmenter}.
 * <p>
 * Every molecule is a combination of structural motifs covering the cases of the fragmentation rules (chains, branched
 * end groups, rings and ring systems, aromatics, conjugated systems, isolated multiple bonds, aromatic remnants as
 * produced by the Ertl fragmenter, stereo centers and stereo double bonds), joined by linkers that are themselves
 * relevant for the rules (direct bond, chains of different lengths, tertiary and quaternary carbons, adjacent branching
 * carbons, multiple bonds). All ordered combinations "motif - linker - motif" and a reproducible random sample of
 * "motif - linker - motif - linker - motif" combinations are generated. The set is generated in code, so no external
 * data file is needed, and the motif names of every molecule tell which cases it combines.
 * </p>
 * <p>
 * The SMILES strings are concatenated, so a motif is bonded to the following part by its last atom and to the
 * preceding part by its first atom. All motifs are chosen so that every atom keeps a valid valence in any position
 * (first, middle or last part); bracket atoms with fixed hydrogen counts are therefore never at a motif's ends. The molecules are parsed like fragments in a MORTAR pipeline: with kekulization or,
 * if that fails (aromatic remnants), without.
 * </p>
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
final class AlkylTestMoleculeGenerator {
    /**
     * Structural motifs, name to SMILES.
     */
    static final Map<String, String> MOTIFS = new LinkedHashMap<>();
    /**
     * Linkers between two motifs, name to SMILES.
     */
    static final Map<String, String> LINKERS = new LinkedHashMap<>();
    /**
     * Number of randomly sampled three-motif combinations.
     */
    static final int TRIPLE_SAMPLE_SIZE = 1500;
    /**
     * Seed of the random sample, fixed for reproducibility.
     */
    static final long SEED = 42;

    static {
        //chains and branched end groups
        AlkylTestMoleculeGenerator.MOTIFS.put("methyl end group", "C");
        AlkylTestMoleculeGenerator.MOTIFS.put("propyl chain", "CCC");
        AlkylTestMoleculeGenerator.MOTIFS.put("octyl chain (longer than default max chain length)", "CCCCCCCC");
        AlkylTestMoleculeGenerator.MOTIFS.put("isopropyl (tertiary C end group)", "C(C)C");
        AlkylTestMoleculeGenerator.MOTIFS.put("tert-butyl (quaternary C end group)", "C(C)(C)C");
        //rings and ring systems
        AlkylTestMoleculeGenerator.MOTIFS.put("cyclopropane ring", "C1CC1");
        AlkylTestMoleculeGenerator.MOTIFS.put("cyclohexane ring", "C1CCCCC1");
        AlkylTestMoleculeGenerator.MOTIFS.put("decalin (fused ring system)", "C1CCC2CCCCC2C1");
        AlkylTestMoleculeGenerator.MOTIFS.put("spiro[4.5]decane (spiro ring system)", "C1CCC2(CC1)CCCC2");
        AlkylTestMoleculeGenerator.MOTIFS.put("norbornane (bridged ring system)", "C1CC2CCC1C2");
        AlkylTestMoleculeGenerator.MOTIFS.put("cyclohexene (ring with double bond)", "C1=CCCCC1");
        AlkylTestMoleculeGenerator.MOTIFS.put("methylenecyclohexane (exocyclic double bond)", "C1CC(=C)CCC1");
        //aromatics
        AlkylTestMoleculeGenerator.MOTIFS.put("benzene ring", "c1ccccc1");
        AlkylTestMoleculeGenerator.MOTIFS.put("naphthalene (fused aromatic)", "c1ccc2ccccc2c1");
        AlkylTestMoleculeGenerator.MOTIFS.put("cyclopentadiene (conjugated ring)", "C1=CC=CC1");
        //conjugated systems
        AlkylTestMoleculeGenerator.MOTIFS.put("butadiene (conjugated D-S-D)", "C=CC=C");
        AlkylTestMoleculeGenerator.MOTIFS.put("hexatriene (conjugated D-S-D-S-D)", "C=CC=CC=C");
        AlkylTestMoleculeGenerator.MOTIFS.put("butadiyne (conjugated T-S-T)", "C#CC#C");
        AlkylTestMoleculeGenerator.MOTIFS.put("butenyne (D-S-T, not conjugated)", "C=CC#C");
        //isolated multiple bonds
        AlkylTestMoleculeGenerator.MOTIFS.put("isolated double bond", "C=C");
        AlkylTestMoleculeGenerator.MOTIFS.put("isolated triple bond", "C#C");
        AlkylTestMoleculeGenerator.MOTIFS.put("allene (cumulated double bonds)", "C=C=C");
        //aromatic remnants (Ertl), the odd ones cannot be kekulized
        AlkylTestMoleculeGenerator.MOTIFS.put("single aromatic remnant atom (no aromatic bond)", "c");
        AlkylTestMoleculeGenerator.MOTIFS.put("aromatic remnant, 2 atoms (kekulizable)", "cc");
        AlkylTestMoleculeGenerator.MOTIFS.put("aromatic remnant, 3 atoms (not kekulizable)", "ccc");
        AlkylTestMoleculeGenerator.MOTIFS.put("aromatic remnant, 5 atoms (not kekulizable)", "ccccc");
        //stereo centers and stereo double bonds
        AlkylTestMoleculeGenerator.MOTIFS.put("tetrahedral stereo center", "C[C@@H](CC)C");
        AlkylTestMoleculeGenerator.MOTIFS.put("E double bond (stereo)", "C/C=C/C");
        AlkylTestMoleculeGenerator.MOTIFS.put("trans-disubstituted cyclohexane (ring stereo centers)", "C[C@H]1CC[C@@H](C)CC1");
        //linkers
        AlkylTestMoleculeGenerator.LINKERS.put("direct bond", "");
        AlkylTestMoleculeGenerator.LINKERS.put("CH2 linker", "C");
        AlkylTestMoleculeGenerator.LINKERS.put("C2 chain linker", "CC");
        AlkylTestMoleculeGenerator.LINKERS.put("C5 chain linker", "CCCCC");
        AlkylTestMoleculeGenerator.LINKERS.put("C9 chain linker (cut by max chain length)", "CCCCCCCCC");
        AlkylTestMoleculeGenerator.LINKERS.put("tertiary C linker", "C(C)");
        AlkylTestMoleculeGenerator.LINKERS.put("quaternary C linker", "C(C)(C)");
        AlkylTestMoleculeGenerator.LINKERS.put("quaternary C linker with ethyl", "C(CC)(C)");
        AlkylTestMoleculeGenerator.LINKERS.put("two adjacent tertiary C", "C(C)C(C)");
        AlkylTestMoleculeGenerator.LINKERS.put("two adjacent quaternary C", "C(C)(C)C(C)(C)");
        AlkylTestMoleculeGenerator.LINKERS.put("double bond linker", "C=C");
        AlkylTestMoleculeGenerator.LINKERS.put("triple bond linker", "C#C");
        AlkylTestMoleculeGenerator.LINKERS.put("C with exo-methylene (double bond branch)", "C(=C)");
        AlkylTestMoleculeGenerator.LINKERS.put("allene linker", "C=C=C");
    }

    /**
     * Private constructor, only static use.
     */
    private AlkylTestMoleculeGenerator() {
    }

    /**
     * Generates the test molecules: all combinations "motif - linker - motif" and a reproducible random sample of
     * "motif - linker - motif - linker - motif". The title of every molecule names its motifs and linkers. SMILES that
     * cannot be parsed are skipped.
     *
     * @return list of the generated molecules
     */
    static List<IAtomContainer> generate() {
        Map<String, String> tmpSmilesByTitle = AlkylTestMoleculeGenerator.generateSmiles();
        List<IAtomContainer> tmpMolecules = new ArrayList<>(tmpSmilesByTitle.size());
        for (Map.Entry<String, String> tmpEntry : tmpSmilesByTitle.entrySet()) {
            IAtomContainer tmpMolecule = AlkylTestMoleculeGenerator.parseAsInPipeline(tmpEntry.getValue());
            if (tmpMolecule != null) {
                tmpMolecule.setTitle(tmpEntry.getKey());
                tmpMolecules.add(tmpMolecule);
            }
        }
        return tmpMolecules;
    }

    /**
     * Generates the SMILES codes of all test molecules (see {@link #generate()}), named by their motifs and linkers.
     *
     * @return map of molecule name to SMILES, in generation order
     */
    static Map<String, String> generateSmiles() {
        List<String> tmpMotifNames = new ArrayList<>(AlkylTestMoleculeGenerator.MOTIFS.keySet());
        List<String> tmpLinkerNames = new ArrayList<>(AlkylTestMoleculeGenerator.LINKERS.keySet());
        Map<String, String> tmpSmilesByTitle = new LinkedHashMap<>();
        for (String tmpFirst : tmpMotifNames) {
            for (String tmpLinker : tmpLinkerNames) {
                for (String tmpSecond : tmpMotifNames) {
                    tmpSmilesByTitle.put(tmpFirst + " + " + tmpLinker + " + " + tmpSecond,
                            AlkylTestMoleculeGenerator.MOTIFS.get(tmpFirst) + AlkylTestMoleculeGenerator.LINKERS.get(tmpLinker) + AlkylTestMoleculeGenerator.MOTIFS.get(tmpSecond));
                }
            }
        }
        Random tmpRandom = new Random(AlkylTestMoleculeGenerator.SEED);
        for (int i = 0; i < AlkylTestMoleculeGenerator.TRIPLE_SAMPLE_SIZE; i++) {
            String tmpFirst = tmpMotifNames.get(tmpRandom.nextInt(tmpMotifNames.size()));
            String tmpFirstLinker = tmpLinkerNames.get(tmpRandom.nextInt(tmpLinkerNames.size()));
            String tmpSecond = tmpMotifNames.get(tmpRandom.nextInt(tmpMotifNames.size()));
            String tmpSecondLinker = tmpLinkerNames.get(tmpRandom.nextInt(tmpLinkerNames.size()));
            String tmpThird = tmpMotifNames.get(tmpRandom.nextInt(tmpMotifNames.size()));
            tmpSmilesByTitle.put(tmpFirst + " + " + tmpFirstLinker + " + " + tmpSecond + " + " + tmpSecondLinker + " + " + tmpThird,
                    AlkylTestMoleculeGenerator.MOTIFS.get(tmpFirst) + AlkylTestMoleculeGenerator.LINKERS.get(tmpFirstLinker) + AlkylTestMoleculeGenerator.MOTIFS.get(tmpSecond)
                            + AlkylTestMoleculeGenerator.LINKERS.get(tmpSecondLinker) + AlkylTestMoleculeGenerator.MOTIFS.get(tmpThird));
        }
        return tmpSmilesByTitle;
    }

    /**
     * Generates a small overview set for manual inspection: every motif bonded by every linker to a methyl end group.
     *
     * @return map of molecule name to SMILES, in generation order
     */
    static Map<String, String> generateOverviewSmiles() {
        Map<String, String> tmpSmilesByTitle = new LinkedHashMap<>();
        for (Map.Entry<String, String> tmpMotif : AlkylTestMoleculeGenerator.MOTIFS.entrySet()) {
            for (Map.Entry<String, String> tmpLinker : AlkylTestMoleculeGenerator.LINKERS.entrySet()) {
                tmpSmilesByTitle.put(tmpMotif.getKey() + " + " + tmpLinker.getKey() + " + methyl end group",
                        tmpMotif.getValue() + tmpLinker.getValue() + "C");
            }
        }
        return tmpSmilesByTitle;
    }

    /**
     * Parses the given SMILES like MORTAR does for fragments in a pipeline: with kekulization and atom type perception
     * or, if that fails, without both.
     *
     * @param aSmiles SMILES to parse
     * @return the parsed molecule, null if the SMILES is invalid
     */
    static IAtomContainer parseAsInPipeline(String aSmiles) {
        try {
            return ChemUtil.parseSmilesToAtomContainer(aSmiles, true, true);
        } catch (CDKException anException) {
            try {
                return ChemUtil.parseSmilesToAtomContainer(aSmiles, false, false);
            } catch (CDKException anInvalidSmilesException) {
                return null;
            }
        }
    }
}
