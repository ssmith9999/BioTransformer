package biotransformer;

import java.io.IOException;
import java.util.ArrayList;

import org.openscience.cdk.CDKConstants;
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
import biotransformer.utils.FileUtilities;

public class HGutBTTest extends HGutBTransformer{

	public HGutBTTest() throws IOException, CDKException{
		// TODO Auto-generated constructor stub
	}

	
	public static void main(String[] args) throws Exception{
		HGutBTTest em = new HGutBTTest();
////		System.out.println(em.getBioSystemName());
		IAtomContainer mol1 = em.getSmiParser().parseSmiles("OC1CC2=C(O)C=C(O)C=C2OC1C3=CC(O)=C(O)C=C3");
		mol1.setProperty(CDKConstants.TITLE, "(+/-)-Epicatechin");
//		IAtomContainer mol2 = em.getSmiParser().parseSmiles("OC1=CC=C(C=C1)C1OC2=CC(O)=CC(O)=C2CC1=O");
//		mol2.setProperty(CDKConstants.TITLE, "5,7‐dihydroxy‐2‐(4‐hydroxyphenyl)‐2,4‐dihydro‐1‐benzopyran‐3‐one");
//		IAtomContainer mol3 = em.getSmiParser().parseSmiles("OC1CC2=C(O)C=C(O)C=C2OC1C1=CC=C(O)C=C1");
//		mol3.setProperty(CDKConstants.TITLE, "3,5,7,4'‐tetrahydroxyflavan");
		IAtomContainerSet mols = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);

		mols.addAtomContainer(mol1);
//		mols.addAtomContainer(mol2);
//		mols.addAtomContainer(mol3);
//		
////		ArrayList<Biotransformation> biotransformations = em.simulateGutMicrobialMetabolism(mol2, true, true, 2, 0.0);
//		ArrayList<Biotransformation> biotransformations = em.applyGutMicrobialMetabolismHydrolysisAndReductionChain(mols, true, true, 8, 0.5);
//		IAtomContainerSet acMetabolites = em.extractAtomContainer(biotransformations);
//		
//		for(IAtomContainer a : acMetabolites.atomContainers()){
//			System.out.println(em.smiGen.create(a));
//		}
//		System.out.println(acMetabolites.getAtomContainerCount());
//		System.out.println(biotransformations.size());
////		for(Biotransformation b : biotransformations){
////			b.display();
////			System.out.println();
////		}
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/test-new-hgut.sdf");
		
		ArrayList<Biotransformation> biotransformations = em.applyGutMicrobialMetabolismHydrolysisAndReductionChain(mols, true, true, 1, 0.5);
		em.saveBioTransformationProductsToSdf(biotransformations, "data/test-new-hgut.sdf");
	
	}	
	


}
