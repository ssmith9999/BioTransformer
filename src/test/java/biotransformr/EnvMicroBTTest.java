package biotransformr;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.ChemStructureManipulator;
import biotransformer.utils.FileUtilities;

public class EnvMicroBTTest extends EnvMicroBTransformer{

	public EnvMicroBTTest() throws IOException, ParseException, CDKException{
		// TODO Auto-generated constructor stub
		super();
	}

	// Ciprofloxacin: c1c2c(cc(c1F)N3CCNCC3)n(cc(c2=O)C(=O)O)C4CC4
	// Diazinon: S=P(OCC)(OCC)Oc1nc(nc(c1)C)C(C)C
	// Chlorpyrifos: CCOP(=S)(OCC)OC1=NC(Cl)=C(Cl)C=C1Cl
	// Imidacloprid : [O-][N+](=O)NC/1=N/CCN\1Cc2cnc(Cl)cc2
	// Disulfoton: S=P(OCC)(SCCSCC)OCC
	// ampicillin : CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)C3=CC=CC=C3)C(=O)N2[C@H]1C(O)=O
	// Nitroglycerin:C(C(CO[N+](=O)[O-])O[N+](=O)[O-])O[N+](=O)[O-]
	public static void main(String[] args) throws Exception{
		EnvMicroBTTest em = new EnvMicroBTTest();
		System.out.println(em.getBioSystemName());
		System.out.println(em.getReactionsList().size());
		
		
		IAtomContainer mol = em.getSmiParser().parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1.Cl");
		IAtomContainer stac = ChemStructureManipulator.standardizeMoleculeWithCopy(mol, true);
		IAtomContainerSet filteredCompounds = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		filteredCompounds.addAtomContainer(mol);
//		System.out.println(em.smiGen.create(stac));
//		ArrayList<Biotransformation> biotransformations = em.applyEnvMicrobialTransformationsChain(stac, true, true, 1);
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
//		em.saveBioTransformationProductsToSdf(biotransformations, "data/BTransformerPaperTest/Nitroglycerin_envMicro_metabolites.sdf");	
//	
//		IAtomContainerSet containers = FileUtils.parseSdf("data/evaluation2/Lambda-cyhalothrin.sdf");
		em.simulateEnvMicrobialDegradationAndSaveToSDF(stac, true, true, 2, 1.0, "/Users/yandj/bla.sdf", true);
		
//		BufferedWriter bw0 = new BufferedWriter(new FileWriter("data/BTransformerPaperTest/Nitroglycerin_envMicro_metabolites.tsv"));		
//		
//		bw0.write("InChI" + "\t" + "InChIKey" + "\t" + "Major Isotope Mass" + "\t" + "\t" + 
//				"Metabolite ID" + "\t" + "\t" + "Precursor InChIKey");
//		bw0.write("\n");
//		for(IAtomContainer a : metabolites.atomContainers()){
//			bw0.write(a.getProperty("InChI") + "\t" + a.getProperty("InChIKey") + "\t" +
//					a.getProperty("Major Isotope Mass") + "\t" + a.getProperty("\t") + "\t" + 
//					a.getProperty("Metabolite ID") + "\t" + a.getProperty("\t") + "\t" + 
//					a.getProperty("Precursor InChIKey"));
//			bw0.write("\n");
//		}
//		bw0.close();
	}
	
	
//	public static void main(String[] args) throws Exception{
//		EnvMicroBTransformer em = new EnvMicroBTransformer();
//		System.out.println(em.getBioSystemName());
//		System.out.println(em.getReactionsList().size());
//		IAtomContainer mol = em.smiParser.parseSmiles("O=C(C)C1=COC=CC1=O");
//		IAtomContainer stac = em.standardizeMolecule(mol, true);
//		System.out.println(em.smiGen.create(stac));
////		IAtomContainerSet acMetabolites = em.applyEnvMicrobialTransformations(stac, true, true);
//		ArrayList<Biotransformation> acMetabolites = em.applyEnvMicrobialTransformationsChain(stac, true, true, 3);
////		System.out.println("Nr. of metabolites: " + acMetabolites.getAtomContainerCount());
////		for(IAtomContainer c : acMetabolites.atomContainers()){
////			System.out.println(em.smiGen.create(c));
////			System.out.println(c.getProperties());
////		}
//	}

}
