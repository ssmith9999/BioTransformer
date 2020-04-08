package biotransformer;

import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.transformation.Biotransformation;
import biotransformer.utils.BiotransformerSequence;
import biotransformer.utils.HumanSuperBioTransformer;
import biotransformer.utils.UniversalBioTransformer;

public class BiotransformerSequenceTest{ 
	
	protected IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
	protected SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
	protected SmilesParser	smiParser		= new SmilesParser(builder);
	
	public BiotransformerSequenceTest() {
		
	}
	public static void main(String[] args) throws Exception {
		BiotransformerSequenceTest btst =  new BiotransformerSequenceTest();
		BiotransformerSequence btseq = new BiotransformerSequence("hgut:1; phaseII:1");
		IAtomContainer ac = btst.smiParser.parseSmiles("COC1=CC(\\C=C\\C(=O)CC(=O)\\C=C\\C2=CC=C(O)C(OC)=C2)=CC=C1O");
//		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		ArrayList<Biotransformation> biots = btseq.runSequence(ac, 0.0);
		UniversalBioTransformer ubt = new UniversalBioTransformer();
		ubt.saveBioTransformationProductsToSdf(biots, "../sequence-test.sdf", false);
	}
	
}
