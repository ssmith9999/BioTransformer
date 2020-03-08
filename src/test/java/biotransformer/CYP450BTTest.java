 package biotransformer;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biomolecule.Enzyme;
import biotransformer.biomolecule.Enzyme.EnzymeName;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.fingerprint.ChemStructureFingerprinter;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.ChemicalClassFinder;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.Utilities;

public class CYP450BTTest extends Cyp450BTransformer{

	public CYP450BTTest(BioSystemName bioSName) throws JsonParseException, JsonMappingException, IOException, CDKException {
		// TODO Auto-generated constructor stub
		super(bioSName);
	}

	public static void main(String[] args) throws Exception {
		CYP450BTTest hCyp450 = new CYP450BTTest(BioSystemName.HUMAN);
		IAtomContainer ac = hCyp450.getSmiParser().parseSmiles("CC(=O)N[C@H]1CCC2=CC(=C(C(=C2C3=CC=C(C(=O)C=C13)OC)OC)OC)OC");
//		ArrayList<Biotransformation> cyp450mets= hCyp450.predictCyp450Biotransformations(ac, true, true, 0.5);
		System.out.println("isOligoOrPolysaccharide: " + ChemicalClassFinder.isOligoOrPolysaccharide(ac));

		ArrayList<Biotransformation> cyp450mets= hCyp450.predictCyp450BiotransformationChain(ac, true, true, 1, 0.5);
		System.out.println(cyp450mets.size());
		
//		hCyp450.saveBioTransformationProductsToCSV(cyp450mets, "data/test_cyp450.csv", true);
		hCyp450.saveBioTransformationProductsToSdf(cyp450mets, "/home/yandj/new_test_cyp450.sdf", false);


//		IAtomContainer acp = ChemStructureManipulator.preprocessContainer(ac);
//
//		ArrayList<Biotransformation> bts = hCyp450.metabolizeWithEnzyme(acp, EnzymeName.CYP1A2, false, false, 0.0);
//
//		IAtomContainerSet acMetabolites = hCyp450.extractProductsFromBiotransformations(bts);
//		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println(hCyp450.smiGen.isomeric().create(c));
//		}
	}
	

}
