package biotransformr;

import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import biotransformer.biomolecule.Enzyme.EnzymeName;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.transformation.MetabolicPathway.MPathwayName;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.ChemicalClassFinder;

public class ECBasedBTTest extends ECBasedBTransformer{

	public ECBasedBTTest(BioSystemName bioSName) throws IOException, ParseException, CDKException, URISyntaxException{
		// TODO Auto-generated constructor stub
		super(bioSName);
	}
	
	public static void main(String[] args) throws Exception{
		ECBasedBTTest em = new ECBasedBTTest(BioSystemName.HUMAN);
		System.out.println(em.getBioSystemName());
		System.out.println(em.getReactionsList().get("ecBasedDeconjugations").size());
		System.out.println(em.getReactionsList().get("ecBasedReactionsNonDeconjugative").size());
		
		IAtomContainer mol = em.getSmiParser().parseSmiles("CCCCCCCCCCCCCCCCCC(O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@@H](CO)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCC=CCCCCCC");
		IAtomContainer stmol = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);

		ArrayList<Biotransformation> biotransformations = em.simulateECBasedMetabolismChain(stmol, true, true, 4, 1.0);
		System.out.println("Nr. of biotransformations: " + biotransformations.size());
		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
		System.out.println("Nr. of metabolites: " + acMetabolites.getAtomContainerCount());
		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println();
			System.out.println(em.smiGen.create(c));
			System.out.println(c.getProperties());
		}

		
	}

}
