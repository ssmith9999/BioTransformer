package biotransformer;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.transformation.Biotransformation;
import biotransformer.utils.BiotransformerSequence;
import biotransformer.utils.ChemStructureExplorer;

public class BiotransformerSequenceTest {
	protected IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
	protected SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
	protected SmilesParser	smiParser		= new SmilesParser(builder);
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		
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
	public void test() throws Exception {
		BiotransformerSequence btseq = new BiotransformerSequence("hgut:1; phaseII:1");
		IAtomContainer molecule 	= smiParser.parseSmiles("CCCC(=O)OC");
		ArrayList<Biotransformation> biots = btseq.runSequence(molecule, 0.0);

		assertTrue("There must be at least 1 biotranformation", (biots != null && biots.size()>0));


	}

}
