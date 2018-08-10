package biotransformr;

import java.io.IOException;
import java.util.ArrayList;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.Utilities;
import predicition.P2Filter;

public class PhaseIIBTTest extends Phase2BTransformer{

	public PhaseIIBTTest(BioSystemName bioSName) throws IOException, ParseException, CDKException{
		// TODO Auto-generated constructor stub
		super(bioSName);
	}

		
	public static void main(String[] args) throws Exception{
		PhaseIIBTTest em = new PhaseIIBTTest(BioSystemName.HUMAN);
//		System.out.println(em.getBioSystemName());
//		
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Glucuronidation")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Sulfonation")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}		
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Acetylation")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Glutathione Transfer")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Glycine Transfer")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}
//		for(MetabolicReaction mr : em.reactionsByGroups.get("Methylation")){
//			System.out.println("queries.put(\"" + mr.getReactionName() + "\" , \"" + mr.getReactantSMARTS().get(0) + "\");");
//		}		
		

		
		P2Filter p2f = new P2Filter();		
		IAtomContainer mol = em.getSmiParser().parseSmiles("[H]C1(OC2=CC(O)=C3CC(=O)[C@]([H])(OC3=C2)C2=CC=CC=C2)OC([H])(C(O)=O)C([H])(O)C([H])(O)C1([H])O");
		IAtomContainer stac = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);		
		ArrayList<String> res = p2f.filter(stac);
		
		System.out.println(res.get(0));
		System.out.println("isInvalidPhase2Metabolite: " + ChemStructureExplorer.isInvalidPhase2Metabolite(mol));

//		em.applyReactionChainFromSdfToSingleSdf("data/test/19-oxotestosterone_EC_based_metabolites_PhaseII_metabolites.sdf",
//				true, true, true, 1, 0.0, "data/test/19-oxotestosterone_EC_based_metabolites_PhaseII_metabolites.sdf");
//
//		
//		em.applyReactionChainFromSdfToSingleSdf("data/Fertaric acid_BioT_sim_metabolites.sdf",
//				true, true, true, 1, 0.0, "data/Fertaric_acid_BioT_sim_metabolites_phaseII-Metabolites.sdf", true);
		
		
		ArrayList<Biotransformation> biotransformations = em.applyPhase2TransformationsChainAndReturnBiotransformations(mol, true, true, true, 1, 0.0);
		System.out.println(biotransformations.size());
		System.out.println(Utilities.selectUniqueBiotransformations(biotransformations).size());
//		em.applyReactionChainFromSDF("data/test/a_Phase_II_metabolites.sdf",true, true, 1, 0.0);
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
//		System.out.println("Nr. of metabolites: " + acMetabolites.getAtomContainerCount());
//		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println(em.smiGen.create(c));
////			System.out.println(c.getProperties());
//		}		
		em.saveBioTransformationProductsToSdf(biotransformations, "data/test_PhaseII_metabolites.sdf");	
//		em.saveBioTransformationsToSDF(biotransformations, "data/test/1,3,5-trinitrobenzene_PhaseII_metabolites.sdf");	
		
//		em.applyReactionChainFromSDF("data/thymol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/nootkatone_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/perillyl_alcohol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/pinene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/pulegone_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		
//		em.applyReactionChainFromSDF("data/terpinen-4-ol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/geraniol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);		
//		em.applyReactionChainFromSDF("data/limonene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);			
//		em.applyReactionChainFromSDF("data/linalool_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/menthol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/myrcene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/1,4-cineole_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/1,8-cineole_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/citral_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/citronellal_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/cuminaldehyde_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/fenchone_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/p-cymene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/carvacrol_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/carvone_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/caryophyllene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
//		em.applyReactionChainFromSDF("data/camphene_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);	
//		em.applyReactionChainFromSDF("data/camphor_CYP450_metabolites.sdf",
//				true, true, 1, 0.0);
		
//		System.out.println( "Size: " + em.getReactionsList().get("Methylation").size());
//		IAtomContainer mol = em.smiParser.parseSmiles("OC1CC2=C(OC1C1=CC(O)=C(O)C=C1)C=C(O)C=C2O");
//
//		IAtomContainer stac = em.standardizeMolecule(mol, true);
//		System.err.println("is metabolizable polyphenol: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(stac));
//		System.err.println("Is invalid? " + ChemStructureExplorer.isInvalidPhaseIIMetabolite(stac));
//		System.out.println(em.smiGen.create(stac));
//////		ArrayList<Biotransformation> biotransformations = em.applyPhaseIITransformations(stac, true, true);
//		ArrayList<Biotransformation> biotransformations = em.applyPhaseIITransformationsChainAndReturnBiotransformations(stac, true, true, 2);
//		System.out.println("Nr. of biotransformations: " + biotransformations.size());
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
//		System.out.println("Nr. of metabolites: " + acMetabolites.getAtomContainerCount());
//		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println(em.smiGen.create(c));
//			System.out.println(c.getProperties());
//		}
//		System.out.println("\n");
//		for(Biotransformation b : biotransformations){
//			b.display();
//		}
//		
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/catechin_PhaseII_metabolites.sdf");	
	}

}
