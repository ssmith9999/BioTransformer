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
		
		
//		IAtomContainer atc = em.smiParser.parseSmiles("CCCCCCCCCCCCCCCC(=O)[O-]");
//		
////		IAtomContainer atc = (IAtomContainer) ChemStructureExplorer.createAtomContainerFromSmilesWithProperties("CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCC", true).get("atomContainer");
////		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(atc);
////		AtomContainerManipulator.convertImplicitToExplicitHydrogens(atc);
//		IAtomContainer atc_2 = ChemStructureManipulator.standardizeMoleculeWithCopy(atc);
//		System.out.println("SMILES OF THE GENERATED OBJECT: " + em.smiGen.isomeric().create(atc_2));
//		
//
//		
//		
////		IAtomContainer mol = em.smiParser.parseSmiles("CC(=CCO[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O)C=CC=C(/C)C=CC1=C(C)CCCC1(C)C");
////		IAtomContainer mol = em.smiParser.parseSmiles("COc1ccc(cc1O)-c1cc(=O)c2c(O)cc(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)cc2o1");
////		IAtomContainer mol = em.getSmiParser().parseSmiles("[H][C@](N)(COP(O)(=O)OC[C@@]([H])(COC(=O)CCCC=C/CC=C/CC=C/CC=C/CCCCC)OC(=O)CCCCCCCC=C/CC=C/CCCCC)C(O)=O");
//		IAtomContainer mol = em.getSmiParser().parseSmiles("CCCCCCCCCCCCCCCCCC(O)[C@H](CO[C@@H]1O[C@H](CO)[C@@H](O[C@H]2O[C@@H](CO)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCC=CCCCCCC");
//		IAtomContainer stmol = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stmol);
////		System.out.println("isInvalidPhaseIIMetabolite: " + ChemStructureExplorer.isInvalidPhaseIIMetabolite(stmol));
////		System.out.println("isMetabolizablePolyphenolOrDerivative: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(stmol));
//		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		molecules.addAtomContainer(stmol);
		
		em.simulateECBasedMetabolismAndSaveToSDF("data/test.sdf", 1, 1.0, "data/", true);
//		
		
//		
//		System.out.println(em.smiGen.isomeric().create(stmol));
//		
//		System.out.println("Is Valid Substrate: " + em.esspredictor.isValidSubstrate(stmol, EnzymeName.EC_3_5_2_X));
//		System.out.println("isInvalidPhaseIIMetabolite: " + ChemStructureExplorer.isInvalidPhaseIIMetabolite(stmol));
		
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(stmol);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stmol);
//		
		
//		ArrayList<Biotransformation> biotransformations = em.simulateECBasedMetabolismChain(stmol, true, true, 4, 1.0);
//		ArrayList<Biotransformation> biotransformations = em.simulateECBasedMetabolism(stmol, true, true, 1.0);
//		ArrayList<Biotransformation> biotransformations = em.metabolizeWithEnzyme(mol, em.enzymesByreactionGroups.get("ecBasedDeconjugations").get(1), true, true, 0.0);
		
//		ArrayList<Biotransformation> biotransformations = em.applyEcBasedTransformationsChain(mol, true, true, 1);
//		ArrayList<Biotransformation> biotransformations = em.applyEcBasedTransformations(mol, true, true, 0);
		
//		ArrayList<Biotransformation> biotransformations = em.applyPathwaySpecificBiotransformations(mol,
//				MPathwayName.GLYCEROPHOSPHOLIPID_METABOLISM, true,true,0.0);
		
//		IAtomContainer mol = em.smiParser.parseSmiles("CCCCCCCCCC(=O)OC(COCCCC)COP(O)(O)=O");
//		IAtomContainer stmol = em.standardizeMolecule(mol, true);
//		System.out.println(em.smiGen.create(stmol));
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//		System.out.println(em.smiGen.create(mol));
//		System.out.println("N_DEALKYLATION_OF_PHOSPHORAMIDE: " +  ChemStructureExplorer.compoundMatchesReactionConstraints(em.reactionsHash.get(ReactionName.N_DEALKYLATION_OF_PHOSPHORAMIDE.toString()),stmol));
//		IAtomContainerSet acMetabolites = em.generateAllMetabolitesFromAtomContainer(mol, 
//				em.reactionsHash.get(ReactionName.BETA_D_GALACTOSIDE_GALACTO_HYDROLYSIS_PATTERN1), true);
//		System.out.println(em.enzymesByreactionGroups.get("ecBasedDeconjugations").get(1).getName());
//		System.out.println(em.enzymesByreactionGroups.get("ecBasedDeconjugations").get(1).getReactionSet().get(0).getReactionName());
//		IAtomContainerSet acMetabolites = em.generateAllMetabolitesFromAtomContainer(stmol, 
//				em.enzymesByreactionGroups.get("ecBasedDeconjugations").get(1).getReactionSet().get(0), true);
//		System.out.println(acMetabolites.getAtomContainerCount());
		
//		System.out.println(ChemicalClassFinder.isSphingoLipid(mol));
//		ArrayList<Biotransformation> biotransformations = em.applyPathwaySpecificBiotransformationsChain(mol, MPathwayName.SPHINGOLIPID_METABOLISM, true, true, 6,0.0);
		
//		System.out.println("Nr. of biotransformations: " + biotransformations.size());
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
//		System.out.println("Nr. of metabolites: " + acMetabolites.getAtomContainerCount());
//		for(IAtomContainer c : acMetabolites.atomContainers()){
////			System.out.println();
//			System.out.println(em.smiGen.create(c));
//			System.out.println(c.getProperties());
//		}
////		System.out.println("\n");
////		for(Biotransformation b : biotransformations){
////			b.display();
////		}
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/test/LacCer_EC_Metabolites.sdf");
	
		
//		for(MetabolicReaction m : em.reactionsByGroups.get("ecBasedReactionsNonDeconjugative")){
//			if(ChemStructureExplorer.compoundMatchesReactionConstraints(m, stmol)){
//				System.out.println(m.name);
//			}
//		}
		
	}

}
