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

package de.unijena.cheminf.mortar.model.filter;

import de.unijena.cheminf.mortar.model.io.Importer;
import de.unijena.cheminf.mortar.model.settings.SettingsContainer;

import org.openscience.cdk.AtomContainerSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Test helper that reads molecules like the MORTAR importer does (SMILES via a non-kekulising CDK SMILES parser, as in
 * the SMILES file import, and SD records via the iterating SDF reader) and applies the importer's preprocessing
 * ({@link Importer#preprocessMoleculeSet(IAtomContainerSet, boolean)}: atom type perception, implicit hydrogen atoms,
 * minimal explicit hydrogen atoms, kekulization), so that the molecule filters are tested on the same kind of atom
 * containers they get at import.
 *
 * @author Manuel Schauer
 * @version 1.0.0.0
 */
class PreprocessingTestImporter extends Importer {
    //<editor-fold desc="Static initializer">
    /**
     * Sets the default locale to British English, required for the message texts of the filter enums.
     */
    static {
        Locale.setDefault(Locale.of("en", "GB"));
    }
    //</editor-fold>
    //
    //<editor-fold desc="Constructor">
    /**
     * Constructor with default global settings.
     */
    PreprocessingTestImporter() {
        super(new SettingsContainer());
    }
    //</editor-fold>
    //
    //<editor-fold desc="Methods">
    /**
     * Parses the given SMILES strings without kekulization, names the molecules "mol0", "mol1", ... with the importer's
     * name property, and preprocesses them like the importer with filling of open valences (MORTAR default).
     *
     * @param aSmilesArray SMILES strings
     * @return preprocessed molecules in the given order
     * @throws CDKException if a SMILES string cannot be parsed
     */
    List<IAtomContainer> parseAndPreprocess(String... aSmilesArray) throws CDKException {
        SmilesParser tmpSmilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        tmpSmilesParser.kekulise(false);
        IAtomContainerSet tmpSet = new AtomContainerSet();
        for (int i = 0; i < aSmilesArray.length; i++) {
            IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(aSmilesArray[i]);
            tmpMolecule.setProperty(Importer.MOLECULE_NAME_PROPERTY_KEY, "mol" + i);
            tmpSet.addAtomContainer(tmpMolecule);
        }
        this.preprocessMoleculeSet(tmpSet, true);
        List<IAtomContainer> tmpList = new ArrayList<>(tmpSet.getAtomContainerCount());
        for (IAtomContainer tmpMolecule : tmpSet.atomContainers()) {
            tmpList.add(tmpMolecule);
        }
        return tmpList;
    }

    /**
     * Reads the first record of the given SD file content with the iterating SDF reader as the importer does, sets the
     * given name with the importer's name property, and preprocesses it like the importer.
     *
     * @param anSdfString SD file content
     * @param aName molecule name to set
     * @return preprocessed molecule
     * @throws IOException if the reader cannot be closed
     */
    IAtomContainer readAndPreprocessSdfRecord(String anSdfString, String aName) throws IOException {
        IAtomContainer tmpMolecule;
        try (IteratingSDFReader tmpReader = new IteratingSDFReader(new StringReader(anSdfString),
                SilentChemObjectBuilder.getInstance())) {
            tmpMolecule = tmpReader.next();
        }
        tmpMolecule.setProperty(Importer.MOLECULE_NAME_PROPERTY_KEY, aName);
        return this.preprocess(tmpMolecule);
    }

    /**
     * Preprocesses the given molecule like the importer.
     *
     * @param aMolecule the molecule
     * @return the same, now preprocessed molecule
     */
    IAtomContainer preprocess(IAtomContainer aMolecule) {
        IAtomContainerSet tmpSet = new AtomContainerSet();
        tmpSet.addAtomContainer(aMolecule);
        this.preprocessMoleculeSet(tmpSet, true);
        return aMolecule;
    }
    //</editor-fold>
}
