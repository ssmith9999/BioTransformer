 package biotransformr;

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
import biotransformer.utils.FileUtilities;
import biotransformer.utils.Utilities;

public class CYP450BTTest extends Cyp450BTransformer{

	public CYP450BTTest(BioSystemName bioSName) throws JsonParseException, JsonMappingException, IOException, CDKException {
		// TODO Auto-generated constructor stub
		super(bioSName);
	}

	public static void main(String[] args) throws Exception {
		CYP450BTTest hCyp450 = new CYP450BTTest(BioSystemName.HUMAN);
		// Promazine: CN(C)CCCN1C2=CC=CC=C2SC2=CC=CC=C12
		// Caffeine: CN1C=NC2=C1C(=O)N(C)C(=O)N2C
		// Omeprazole: 
		IAtomContainer ac = hCyp450.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		ArrayList<Biotransformation> cyp450mets= hCyp450.predictCyp450Biotransformations(ac, true, true, 0.5);
		System.out.println(cyp450mets.size());
//		System.out.println("SAME BIOTRANSFORMATION? " + cyp450mets.get(0).equals(cyp450mets.get(0)));
//		System.out.println("SAME BIOTRANSFORMATION? " + cyp450mets.get(0).equals(cyp450mets.get(1)));
//		ac.setProperty(CDKConstants.TITLE,"");
//		IAtomContainer acp = ChemStructureManipulator.preprocessContainer(ac);
//		System.out.println(hCyp450.smiGen.create(acp));
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(acp);
//		System.out.println(hCyp450.esspredictor.isValidCyp450Substrate(a, EnzymeName.CYP2C9));
//		ArrayList<Biotransformation> bts = hCyp450.metabolizeWithEnzyme(acp, EnzymeName.CYP1A2, false, false, 0.0);
//		hCyp450.saveBioTransformationProductsToSdf(cyp450mets, "/Users/yandj/epicatechin-cyp450-mets-1.sdf");
//		IAtomContainerSet containers = hCyp450.extractAtomContainer(bts);
//		String res = FileUtilities.saveAtomContainersToString(containers);
//		System.out.println(res);
//		System.out.println(FileUtilities.saveAtomContainerToString(ac));
		
//		int i = 0;
//		for (MetabolicReaction metR : hCyp450.reactionsHash.values()){
//			System.out.println("queries.put(\"" + metR.name.toLowerCase() + "\", \"" + metR.getReactantSMARTS().get(0) + "\");");
////			System.out.println(metR.name + "\t" + metR.getReactantSMARTS().get(0));
//			i++;
//		}
//		System.out.println(i);
		
//		System.out.println(ChemStructureFingerprinter.generateBTRailsCountFingerprint(ac));
//		System.out.println(ChemStructureExplorer.ringCount(ac));
//		LinkedHashMap<Object, Object> props =  new LinkedHashMap<Object, Object>();
//		props.put("Synonym", "test");
//		Utilities.annotateAtomContainerWithProps(acp, props);
//		System.out.println(acp.getProperty(CDKConstants.TITLE));
//		System.out.println(acp.getProperty("Major Isotope Mass"));
//		
//		LinkedHashMap<String, Object> convers = ChemStructureExplorer.createAtomContainerFromSmilesWithProperties("CCN(CC)CCOC(=O)C1(CCCCC1)C2CCCCC2", true);
//		System.out.println(convers);
		
//		System.out.println(hCyp450.esspredictor.isValidSubstrate(a, EnzymeName.CYP2C9));
		
//		ArrayList<Biotransformation> bts = hCyp450.applyReactionAndReturnBiotransformations(ac,
//				hCyp450.bSystem.getReactionsHash().get("AROMATIC_HYDROXYLATION_OF_FUSED_BENZENE_RING_PATTERN2"), true, 0.0);
//		System.out.println("Macthing? " + ChemStructureExplorer.compoundMatchesReactionConstraints(
//				hCyp450.bSystem.getReactionsHash().get("AROMATIC_HYDROXYLATION_OF_FUSED_BENZENE_RING_PATTERN2"), a));
		
//		ArrayList<Biotransformation> bts = hCyp450.predictCyp450Biotransformations("data/Chlorpyrifos.sdf", true, true, 0.0);		

//		IAtomContainerSet acMetabolites = hCyp450.extractAtomContainer(bts);
//		for(IAtomContainer c : acMetabolites.atomContainers()){
//			System.out.println(hCyp450.smiGen.isomeric().create(c));
//			System.out.println("\t" + c.getProperty("Reactions"));
//		}
//		hCyp450.saveBioTransformationProductsToSdf(bts, "data/BTransformerPaperTest/Chlorpyrifos_cyp450_metabolites.sdf");
		
		
//		LinkedHashMap<ReactionName, Set<EnzymeName>> rtoenz = new LinkedHashMap<ReactionName, Set<EnzymeName>>();
//		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
//		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader("/Users/yandj/Programming/Projects/Metabolism/reactantpredictor/data/JarleiSet_CYP_predictions.sdf"),
//				bldr);
//		
////		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader("/Users/yandj/Programming/Projects/Metabolism/bioT/biotransformr/data/1_8-cineole.sdf"), bldr);
//		
//		while (sdfr.hasNext()){
//			IAtomContainer mol_ = sdfr.next();
//			IAtomContainer mol = ChemStructureManipulator.preprocessContainer(mol_);
////			System.out.println("New Molecule: " + mol.getProperty(CDKConstants.TITLE));
//			ArrayList<Biotransformation> bts = new ArrayList<Biotransformation>();
////			try{
//				
////			String[] enzlist;
////			String[] reacList
////			System.out.println(mol.getProperty("1A2").toString().contentEquals("R"));
//			
//			if(mol.getProperty("1A2").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP1A2");
////				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP1A2, false, false, 0.0));
////				}catch(NullPointerException n){
////					System.err.println(n.getMessage());
////				}
//			}
////			System.out.println(bts.size());
//			
//			if(mol.getProperty("2A6").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2A6");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2A6, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			if(mol.getProperty("2B6").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2B6");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2B6, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2C8");
//			if(mol.getProperty("2C8").toString().contentEquals("R")){
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2C8, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			if(mol.getProperty("2C9").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2C9");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2C9, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			if(mol.getProperty("2C19").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2C19");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2C19, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			if(mol.getProperty("2D6").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2D6");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2D6, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			
//			if(mol.getProperty("2E1").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP2E1");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP2E1, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//			
//			
//			if(mol.getProperty("3A4").toString().contentEquals("R")){
//				System.out.println(mol.getProperty(CDKConstants.TITLE) + " <==> CYP3A4");
//				try{
//				bts.addAll(hCyp450.metabolizeWithCyp450(mol, EnzymeName.CYP3A4, false, false, 0.0));
//				}catch(NullPointerException n){
//					System.err.println(n.getMessage());
//				}
//			}
//				
//			ArrayList<Biotransformation> cyp450 = hCyp450.predictCyp450Biotransformations(ac, true, true, 0.0);			
//			System.out.println("bts.size() " + cyp450.size());
//			if(cyp450.size()>0){
//				for(Biotransformation bt : cyp450){
//					System.out.println(hCyp450.smiGen.create(bt.getProducts().getAtomContainer(0)));
//				}
//				hCyp450.saveBioTransformationProductsToSdf(cyp450, "data/test/" + ac.getProperty(CDKConstants.TITLE).toString().replaceAll(" ", "_") + "_CYP450_metabolites.sdf");			
//			}
////			}
////			catch(Exception e) {
////				System.err.println("ERROR");
////				System.err.println(e.getMessage());
////			}
//			
//			}
		
//			IAtomContainerSet metabolites = FileUtils.parseSdf("data/test/CID22096324-test.sdf");
//			hCyp450.simulateCyp450MetabolismAndSaveToSDF(metabolites, 1, 1.0, "data/test/");
			
//			BufferedWriter bw0 = new BufferedWriter(new FileWriter("data/BTransformerPaperTest/Disulfoton_CYP450_metabolites.tsv"));		
//			
//			bw0.write("InChI" + "\t" + "InChIKey" + "\t" + "Major Isotope Mass" + "\t" + "\t" + 
//					"Metabolite ID" + "\t" + "\t" + "Precursor InChIKey");
//			bw0.write("\n");
//			for(IAtomContainer a : metabolites.atomContainers()){
//				bw0.write(a.getProperty("InChI") + "\t" + a.getProperty("InChIKey") + "\t" +
//						a.getProperty("Major Isotope Mass") + "\t" + a.getProperty("\t") + "\t" + 
//						a.getProperty("Metabolite ID") + "\t" + a.getProperty("\t") + "\t" + 
//						a.getProperty("Precursor InChIKey"));
//				bw0.write("\n");
//			}
//			bw0.close();
			
			
//		ArrayList<Biotransformation> biotransformations = hCyp450.predictCyp450Biotransformations("data/JarleiSet.sdf", true, false, 0.0);
//		hCyp450.saveBioTransformationProductsToSdf(biotransformations,"data/JarleiSet_CYP450_metabolites.sdf");
	
		
//		System.out.println(hCyp450.reactionsByGroups.get("cypReactions").size());
//		int count = 0;
//		for(Enzyme e: hCyp450.enzymesList){
//			System.out.println(e.getName() + " : " + e.getReactionsNames().size());
//			count = count + e.getReactionsNames().size();
//		}
//		System.out.println(count);
	}
	

}
