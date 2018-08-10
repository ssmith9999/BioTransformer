package biotransformr;

import java.io.IOException;
import java.util.ArrayList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSReaction;
import biotransformer.biomolecule.Enzyme;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MReactionSets;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.utils.ChemStructureExplorer;

public class HGutBTTest extends HGutBTransformer{

	public HGutBTTest() throws IOException, CDKException{
		// TODO Auto-generated constructor stub
	}

	
	public static void main(String[] args) throws Exception{
		HGutBTTest em = new HGutBTTest();
		System.out.println(em.getBioSystemName());
//		System.out.println(em.getReactionsList().get("gutMicroReductiveReactions").size());
//		
		em.simulateGutMicrobialMetabolismAndSave("data/test.sdf", true, true, 1, 0.0, "data/test_results.sdf");
		
		// OC(=O)C=CCC1=CC(O)=C(O)C(O)=C1
		// OC(=O)C1=CC(O)=CC2=C(C=C(O)C=C12)C(O)=O
		// OC(=O)C1=CC(O)=C2OC(=O)C3=C2C1=CC(O)=C3
		// C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](OC3=C(OC4=CC(O)=CC(OC(C)=O)=C4C3=O)C3=CC(O)=C(OP(O)([O-])=O)C=C3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@H](O)[C@H]1O
		// catechin: O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC(O)=C(O)C=C1
		// O[C@@H]1Cc2c(O)cc3O[C@@]4(Oc5cc(O)cc(O)c5[C@@H]([C@H]4O)c3c2O[C@@H]1c1ccc(O)c(O)c1)c1ccc(O)c(O)c1
		// Ellagic acid : C1=C2C3=C(C(=C1O)O)OC(=O)C4=CC(=C(C(=C43)OC2=O)O)O
		// Epigallocatechin: C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C(=C3)O)O)O)O
		// Hesteretin: [H]OC1=C([H])C(O[H])=C2C(O[C@]([H])(C3=C([H])C(O[H])=C(OC([H])([H])[H])C([H])=C3[H])C([H])([H])C2=O)=C1[H]
		// Ecpicatechin : O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1
		// OC[C@H]1O[C@@H](OC2=C([O+]=C3C=C(O)C=C(O)C3=C2)C2=CC=C(O)C(O)=C2)[C@H](O)[C@@H](O)[C@@H]1O
		// OC[C@H]1O[C@@H](OC2=CC3=C(C=C(O)C=C3O)[O+]=C2C2=CC=C(O)C(O)=C2)[C@H](O)[C@@H](O)[C@@H]1O
//		IAtomContainer mol = em.getSmiParser().parseSmiles("COC1=C(O)C=CC(=C1)C1=[O+]C2=CC(O)=CC(O)=C2C=C1O");
		IAtomContainer mol2 = em.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
//		System.out.println(ChemStructureExplorer.isBioTransformerValid(mol2));
		
//		ArrayList<Biotransformation> bts = em.applyGutMicrobialReductions(mol2, true, true, 0.0);
//		
//		IAtomContainerSet acs = em.extractAtomContainer(bts);
//		for(IAtomContainer a : acs.atomContainers()){
//			System.out.println(em.smiGen.create(a));
//		}
		
//		em.simulateGutMicrobialPolyphenolMetabolismAndSave("data/Quercetin.sdf", true, true, 1, 1.0, "data/Quercetin-hgut.sdf", false);
		
//		IAtomContainer mol = em.smiParser.parseSmiles("OC[C@H]1OC(OC2=CC(OS([O-])(=O)=O)=C([N-][N+]#N)C(OP([O-])([O-])=O)=C2)[C@H](O)[C@@H](O)[C@@H]1O");				
//		SMIRKSReaction sr = em.smrkMan.parse("[#8;X1-:3][P;X4:2]([#8;X1-:4])(=[O;X1:5])[#8;X2:1]-[*,#1:6]>>[H][#8;X2:3][P;X4:2](=[O;X1:5])([#8;X2:4][H])[#8;X2:1]-[*,#1:6]");
//		
////		System.out.println("Before adding hydrogens: "  + em.smiGen.isomeric().create(mol));
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
////		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//		System.out.println("After adding hydrogens: "  + em.smiGen.isomeric().create(mol));
//		mol = em.smrkMan.applyTransformationWithSingleCopyForEachPos(mol, null, 
//				sr).getAtomContainer(0);
//		System.out.println(em.smiGen.isomeric().create(mol));
		
//		IAtomContainer smol = AtomContainerManipulator.suppressHydrogens(mol);
//		System.out.println("After copying and suppressing hydrogens: "  + em.smiGen.isomeric().create(mol));
//		
//		System.out.println("\n" + MReactionSets.standardizationReactions.get(12).getReactionName());
//		mol = em.smrkMan.applyTransformationWithSingleCopyForEachPos(mol, null, 
//				MReactionSets.standardizationReactions.get(12).getSmirksReaction()).getAtomContainer(0);
//		System.out.println(em.smiGen.isomeric().create(mol));
//		
//		System.out.println("\nMolecules with suppressed hydrogens: "  + em.smiGen.isomeric().create(smol));
//		System.out.println(MReactionSets.standardizationReactions.get(11).getReactionName());
//		smol = em.smrkMan.applyTransformationWithSingleCopyForEachPos(smol, null, 
//				MReactionSets.standardizationReactions.get(11).getSmirksReaction()).getAtomContainer(0);
//		System.out.println(em.smiGen.isomeric().create(smol));
		
		
//		System.out.println("\n" + MReactionSets.standardizationReactions.get(11).getReactionName());
//		
//		
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//		IAtomContainerSet sets = em.smrkMan.applyTransformationWithSingleCopyForEachPos(mol, null, 
//				MReactionSets.standardizationReactions.get(11).getSmirksReaction());
//		System.out.print("sets: ");
//		System.out.println(sets.getAtomContainerCount());
//		mol 
//		System.out.println(em.smiGen.isomeric().create(mol));		
		
		
//		IAtomContainer stac = em.standardizeMoleculeWithCopy(mol, true);
////		
////		AtomContainerManipulator.convertImplicitToExplicitHydrogens(stac);
////		//		IAtomContainer stac = mol;
//		System.out.println("stac: " + em.smiGen.isomeric().create(stac));
//		System.out.println("mol: " +  em.smiGen.isomeric().create(mol));
////
//		
////		IAtomContainer mol2 = em.smiParser.parseSmiles("[H]OC1=C(O[H])C([H])=C(C([H])=C([H])[H])C([H])=C1[H]");
////		System.out.println("MOL2: " + em.smiGen.isomeric().create(mol2));
//		
//		System.out.println("Is metabolizable polyphenol derivative: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(stac));
//		System.out.println("isPolyphenolOrDerivative: " + ChemStructureExplorer.isPolyphenolOrDerivative(stac));
//		System.out.println("isInvalidPhaseIIMetabolite: " + ChemStructureExplorer.isInvalidPhaseIIMetabolite(stac));
		
//		ArrayList<Biotransformation> biotransformations = em.applyGutMicrobialDeconjugationsChain(stac, true, true , 4, 0.0);
//		ArrayList<Biotransformation> biotransformations = em.applyGutMicrobialReductionsChain(stac, true, false , 6, 0.0);

	//		ArrayList<Biotransformation> biotransformations = em.applyReactionsChainAndReturnBiotransformations(stac, em.reactionsList.get("gutMicroReactions"), true, true, 2, 0.0);

//		System.out.println(MReactionSets.standardizationReactions.get(11).getReactantSMARTS().get(0));
//		System.out.println(em.bSystem.getReactionsHash().get("HYDROLYSIS_OF_CARBOXYLIC_ACID_ESTER_PATTERN2").getReactionSMIRKS());
//		ArrayList<Biotransformation> biotransformations = em.applyReactionAtOnceAndReturnBiotransformations(stac,
//				em.bSystem.getReactionsHash().get("DECARBOXYLATION_OF_FUSED_BENZENE"), false, 0.0);
		
//		System.out.println("Match? " + ChemStructureExplorer.compoundMatchesReactionConstraints(
//				em.bSystem.getReactionsHash().get("LACTONE_HYDROLYSIS_PATTERN1"), mol));
//		ArrayList<Biotransformation> biotransformations = em.applyReactionAndReturnBiotransformations(mol,
//				em.bSystem.getReactionsHash().get("LACTONE_HYDROLYSIS_PATTERN1"), true, 0.0);
//		
//		ArrayList<Biotransformation> biotransformations = em.applyReactionAtOnceAndReturnBiotransformations(mol,
//				em.reactionsList.get("deconjugationReactions"), true, true ,4, 0.0);
//		
//		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//		
////		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol2);
////		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol2);
//		
//		IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
////		IAtomContainerSet molecules2 = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
//		molecules.addAtomContainer(mol);
////		molecules.addAtomContainer(mol2);
//		
////		molecules2 = (IAtomContainerSet) molecules.clone();
//		
//		System.err.println("Metabolizable polyphenol: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(mol));
////		ArrayList<Biotransformation> biotransformations = em.simulateGutMicrobialPolyphenolMetabolism(molecules, true, true, 2, 0.0);
//		ArrayList<Biotransformation> biotransformations = em.metabolizeWithEnzymesBreadthFirst(molecules, em.enzymesByreactionGroups.get("gutMicroReductiveReactions"), true, true, 6, 0.0);
////		ArrayList<Biotransformation> biotransformations = em.metabolizeWithEnzymes(molecules, em.enzymesByreactionGroups.get("gutMicroReductiveReactions"), true, true, 0.0);
//		System.out.println("There are a total of biotransformations " + biotransformations.size());
////		ArrayList<Biotransformation> biotransformations2 = em.applyReactionsChainAndReturnBiotransformations(molecules2, em.reactionsByGroups.get("gutMicroReductiveReactions"), true, true, 2, 0.0);
////		System.out.println("There are a total of biotransformations (2) " + biotransformations2.size());
//		
//		
//		int nrReact1=0;
//		for(Enzyme e : em.enzymesByreactionGroups.get("gutMicroReductiveReactions")){
//			nrReact1 = nrReact1 + e.getReactionSet().size();
//		}
//		int nrReact2=em.reactionsByGroups.get("gutMicroReductiveReactions").size();
//		
//		System.out.println("nrReact1 : " + nrReact1);
//		System.out.println("nrReact2 : " + nrReact2);
//		
//		ArrayList<Biotransformation> biotransformations = em.simulateGutMicrobialMetabolism(mol2,
//				true, true, 6, 0.0);
//		
//		System.err.println("IS DECONJUGATION CANDIDATE: " + em.isDeconjugationCandidate(mol2));
		
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/epicatechin-hgut-6.sdf", true);
		
//		ArrayList<Biotransformation> biotransformations = em.applyGutMicrobialReductionsChain(mol, true, true , 8, 0.0);
		
//		MetabolicReaction mreact = em.smrkMan.parse(">>");
		
//		ArrayList<Biotransformation> biotransformations = em.applyReactionAndReturnBiotransformations(mol, em.reactionsHash.
//				get(ReactionName.FLAVANON_3_OL_C_RING_FISSION.toString()), false);
//		System.out.println(biotransformations.size());	
//		IAtomContainerSet r = em.generateAllMetabolitesFromAtomContainer(mol, 
//				em.reactionsHash.get(ReactionName.FLAVONOID_C_RING_REDUCTION.toString()), true);
//		System.out.println(r.getAtomContainerCount());

////////		System.out.println(biotransformations.get(0).getClass());
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
////		
////		for(Biotransformation b : biotransformations){
////			b.display();
////			System.out.println();
////		}
////		
////		System.out.println(acMetabolites.getAtomContainerCount()  + " metabolites from " + biotransformations.size() + " transformations.");
////		
//		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println(em.smiGen.isomeric().create(c));
//			System.out.println("\t" + c.getProperty("Reactions"));
//		}		
//		
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/HGutM_Metabolites_6Steps_1.sdf");	
//		em.saveBioTransformationProductsToSdf(biotransformations2, "data/HGutM_Metabolites_6Steps_2.sdf");
//		System.out.println("There are a total of " + acMetabolites.getAtomContainerCount() + " metabolites");
//		System.out.println("ecBasedReactionsNonDeconjugative: " + em.enzymesByreactionGroups.get("gutMicroPhaseIIReactions").size());
//		System.out.println("ecBasedReactionsNonDeconjugative: " + em.enzymesByreactionGroups.get("gutMicroReductiveReactions").size());

//		System.out.println("Is invalid phase II metabolite: " + ChemStructureExplorer.isInvalidPhaseIIMetabolite(mol));
//		IAtomContainer m = em.smiParser.parseSmiles("[H]OC(=O)C1([H])OC([H])(OC2=C(O[H])C([H])=C(C([H])=C2[H])C2([H])OC3=C(C(O[H])=C([H])C(O[H])=C3[H])C([H])([H])C2([H])O[H])C([H])(O[H])C([H])(O[H])C1([H])O[H]");
//		System.out.println("Is molecule a valid phase II metabolite? " + (!ChemStructureExplorer.isInvalidPhaseIIMetabolite(m)));
//		System.out.println("Is metabolizable polyphenol derivative: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(m));
//	
//		IAtomContainer n = em.smiParser.parseSmiles("C1(=O)C2=CC=CC=C2C3=CC=CC=C3O1");
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(n);
//		System.out.println("Is molecule a valid phase II metabolite? " + (!ChemStructureExplorer.isInvalidPhaseIIMetabolite(n)));
//		System.out.println("Is metabolizable polyphenol derivative: " + ChemStructureExplorer.isMetabolizablePolyphenolOrDerivative(n));
//	
//		
//		SmartsPatternCDK isoflavonePattern = new SmartsPatternCDK("[R0;*,#1]-[#6;R1]=,:1[#8]c2c([#6;R1](=[O;X1])[#6;R1]=,:1-[c;R1]1[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1](-[R0;*,#1])[c;R1]1-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c(-[R0;*,#1])c2-[R0;*,#1]");
//		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");
//		System.out.println("Is flavonoid with 2 or more sulfate groups " + (isoflavonePattern.hasSMARTSPattern(n) == 0 && sulfatedRadicalPattern.hasSMARTSPattern(n) >= 2));

		
//		IAtomContainer ac = em.smiParser.parseSmiles("C1C(C(OC(=O)C2=CC(=C(C(=C2C3=C(C(=C4C5=C3C(=O)OC6=C(C(=C(C7=C(C(=C(C=C7C(=O)O1)O)O)O)C(=C56)C(=O)O4)O)O)O)O)O)O)O)C8C(OC(=O)C9=CC(=C(C(=C9C1=C(C(=C(C=C1C(=O)O8)O)O)O)O)O)O)C=O)O");
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(ac);
//		System.out.println(em.smiGen.create(ac));
//		List<List<Integer>> l = ChemStructureExplorer.findAllOccurences("[H][#8]-[#6;R1]=,:1[#6;R1](-[#8])=,:[#6;R1](-[#8])[#6;R1]=,:[#6]-2[#6]=,:1-[#6]1=,:[#6;R1](-[#8][H])[#6;R1](-[#8])=,:[#6;R1](-[#8])[#6;R1]=,:[#6]1-[#6](=O)-[#8]-[#6]-[#6]-[#8]-[#6]-2=O"
//				,ac);
//		System.out.print("Number of occurrences: " + l.size() + "\n\n\n\n");
////		List<List<Integer>> m = ChemStructureExplorer.findAllOccurences(""
//				,ac);
//		System.out.print(m.size());
		
//		System.out.println( compoundMatchesReactionConstraints(em.reactionsList.get("gutMicroReactions").,
//				ac));
		
	}	
	


}
