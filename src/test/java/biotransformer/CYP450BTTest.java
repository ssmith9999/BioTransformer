package biotransformer;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.ArrayList;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureExplorer;

public class CYP450BTTest {
	static Cyp450BTransformer hCyp450;
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		hCyp450 = new Cyp450BTransformer(BioSystemName.HUMAN);
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test() throws Exception{
		IAtomContainer molecule = hCyp450.getSmiParser().parseSmiles("CC(=O)Nc1ccccc1");
		IAtomContainer acetaminophen = hCyp450.getSmiParser().parseSmiles("CC(=O)NC1=CC=C(C=C1)O");
		
		ArrayList<Biotransformation> cyp450bios = hCyp450.predictCyp450BiotransformationChain(
				molecule, true, true, 1, 0.5);
		IAtomContainerSet cyp450mets = hCyp450.extractProductsFromBiotransformations(cyp450bios);
		
		assertTrue("There must be at least 1 biotranformation", (cyp450bios != null && cyp450bios.size()>0));
		
		assertTrue("Acetaminophen must be a metabolite (CC(=O)NC1=CC=C(C=C1)O)", 
				ChemStructureExplorer.atomContainerInclusionHolds(cyp450mets, acetaminophen));
			
	}

}
