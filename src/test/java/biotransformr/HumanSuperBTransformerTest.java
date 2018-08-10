package biotransformr;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.Time;
import java.util.ArrayList;

import org.json.JSONObject;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.mashape.unirest.http.HttpResponse;
import com.mashape.unirest.http.JsonNode;
import com.mashape.unirest.http.Unirest;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.transformation.MetabolicReaction;
import biotransformer.transformation.MRPatterns.ReactionName;
import biotransformer.utils.ChemStructureExplorer;
import biotransformer.utils.ChemdbRest;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HumanSuperBioTransformer;


public class HumanSuperBTransformerTest extends HumanSuperBioTransformer{

	public HumanSuperBTransformerTest() throws IOException, ParseException, CDKException {
		// TODO Auto-generated constructor stub
		super();
	}
	
//	public static void main(String[] args) throws CDKException {
//	    String             smi = "C[C@]12CC[C@H]3[C@@H](CCC4=CC(O)=CC=C34)[C@@H]1CC[C@@]2(O)C#C";
//	    IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
//	    SmilesParser       smipar = new SmilesParser(bldr);
//	    IAtomContainer mol = smipar.parseSmiles(smi);
//	    System.out.println(InChIGeneratorFactory.getInstance().getInChIGenerator(mol).getInchiKey());
//	    // not needed
//	    // AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//	    AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//	    System.out.println(InChIGeneratorFactory.getInstance().getInChIGenerator(mol).getInchiKey());
//	}

	
	public static void main(String[] args) throws Exception{
		HumanSuperBTransformerTest hsbt 	= new HumanSuperBTransformerTest();
		IChemObjectBuilder 	builder 	= SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		
		// Metformin: CN(C)C(=N)NC(N)=N
		// carvone: 
//	IAtomContainer molecule = smiParser.parseSmiles("OC(=O)CCCCCCCCCCCCCCCCC(O)=O");
//		Biotransformer b = new Biotransformer(BioSystemName.HUMAN);
//		IAtomContainer mt = b.standardizeMoleculeWithCopy(molecule);
//				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mt);
//		ArrayList<Biotransformation> bt = hsbt.simulateHumanMetabolism(mt);
//		
		long t0 = System.currentTimeMillis();
		
		
//		ArrayList.
		
		hsbt.simulateHumanAndGutMicrobialMetabolismFromSDF("/Users/yandj/Quercetin.sdf", true);
		
//		hsbt.simulateHumanMetabolismFromSDFtoSingleSDF("data/INRA/camphor.sdf", "data/INRA/camphor_superbioMetabolites.sdf", true);

//		hsbt.simulateHumanMetabolismFromSDFtoSingleSDF("/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_random_1000.sdf","/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_random_1000_superbio_mets.sdf", false);
		
//		hsbt.simulateHumanMetabolismFromSDFtoSingleSDF("data/Monocaffeoyl(-)-tartaric acid.sdf","data/Monocaffeoyl(-)-tartaric_acid_mets.sdf", false);
		
//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSingleSDF("/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_random_1000.sdf", "/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_random_annotated_1000_metabolites.sdf",1,0.0, true);

//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSingleSDF("/Users/yandj/epicatechin.sdf", "/Users/yandj/epicatechin_AllHuman_1Steps_Metabolites.sdf",1,0.0, true);

//		IAtomContainerSet s = FileUtilities.parseSdf("data/evaluation5/biotransformerEvaluationSet_PharmaceuticalDrugsAndPesticides.sdf");
//		ArrayList<Biotransformation> bts = hsbt.simulateOneStepAllHuman(s.getAtomContainer(0));
		
//		System.out.println("Substrate INCHIKEY" + bts.get(0).getSubstrates().getAtomContainer(0).getProperty("InChIKey"));
		//		
//		long t1 = System.currentTimeMillis();
//		
//		System.out.println("It took: " + (t1 - t0) + " milliseconds to compute.");
				
//		System.out.println("Nr. of unique compounds: " + FileUtils.countUniqueCompounds("/Users/yandj/Projects/DB_data/FooDB/foodb_2017_06_29_csv/compounds_random_annotated_1000_metabolites.sdf"));
		
		
//		hsbt.simulateHumanMetabolismFromSDF("data/Fertaric_acid.sdf","data/",true);
//		hsbt.simulateHumanMetabolismFromSDFtoSingleSDF("data/Carbaryl.sdf","data/Carbaryl.sdf_bioT_supertransformer_metabolites.sdf", true);
//		ArrayList<Biotransformation> bt = hsbt.simulateOneStep(molecule);
//		hsbt.ecb.saveBioTransformationProductsToSdf(bt, "data/test_bioT_metabolites.sdf");
//		hsbt.simulateHumanMetabolismOneStepFromSDFtoSingleSDF("data/evaluation2/17-Ethinylestradiol.sdf", "data/evaluation2/17-Ethinylestradiol_biot_Metabolites_OneStep.sdf");
//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSingleSDF("data/1‐methyl‐6‐(pyridin‐3‐yl)piperidin‐2‐ol.sdf", "data/1‐methyl‐6‐(pyridin‐3‐yl)piperidin‐2‐ol_biot_Metabolites_OneStep.sdf",1, 1.0);
//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSDF("data/Metronidazole.sdf", "data/",1, 1.0, true);
		
//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSDF("data/evaluation5/biotransformerEvaluationSet_PharmaceuticalDrugsAndPesticides.sdf", "data/evaluation5/",1, 1.0, true);
//		hsbt.predictMetabolismAllHumanFromSDFAndSavetoSDF("data/evaluation6/biotransformerEvaluationSet_TestSet2.sdf", "data/evaluation6/",1, 1.0, true);
		
		
//		HttpResponse<JsonNode> jsonResponse = Unirest.post("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/IHWFPRKZRRGTTI-UHFFFAOYSA-N/cids/json")
//				.header("accept", "application/json").asJson();
//		
//		JSONObject jObject = jsonResponse.getBody().getObject();
//		//https://www.programcreek.com/java-api-examples/index.php?api=com.mashape.unirest.http.JsonNode
//		System.out.println(jObject.get("IdentifierList")); // {"CID":[14432748]}
//		
//		JSONObject identifierList = new JSONObject(jObject.get("IdentifierList").toString());
//		
//		System.out.println(identifierList.getJSONArray("CID").get(0));
//	
//		
//		HttpResponse<JsonNode> jsonResponseSyn = Unirest.post("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/6429037/synonyms/json")
//				.header("accept", "application/json").asJson();
//
//		JSONObject jObjectSyn = jsonResponseSyn.getBody().getObject();
//		System.out.println(jObjectSyn.get("InformationList"));
//		
//		JSONObject informationList = new JSONObject(jObjectSyn.get("InformationList").toString());
//		
//		System.out.println(informationList.getJSONArray("Information").getJSONObject(0).getJSONArray("Synonym").get(0));
		
//		HttpResponse<JsonNode> jsonResponse = Unirest.post("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/GRWPTSXPZYCYOM-UHFFFAOYSA-N/synonyms/json")
//				.header("accept", "application/json").asJson();
//		
//		JSONObject jObject = jsonResponse.getBody().getObject();
////		System.out.println(jObject.get("Fault").toString());
//		//https://www.programcreek.com/java-api-examples/index.php?api=com.mashape.unirest.http.JsonNode
//		JSONObject information = new JSONObject(jObject.get("InformationList").toString());	
//		System.out.println(information.getJSONArray("Information").get(0));
//		JSONObject synonyms = new JSONObject(information.getJSONArray("Information").get(0).toString());
//		System.out.println(synonyms.get("CID").toString()); // 2
////		System.out.println(synonyms.getJSONArray("Synonym").get(0));
//		
////		for(int i = 0; i < synonyms.getJSONArray("Synonym").length(); i++){
////			System.out.println(synonyms.getJSONArray("Synonym").get(i));
////		}
		
//		System.out.println(ChemdbRest.getSynonymsObjectViaInChIKey("DYBBEZHAELJFKW-ZLJBKAMKSA-N").get("CID").get(0));
//		System.out.println(ChemdbRest.getSynonymsObjectViaInChIKey("GRWPTSXPZYCYOM-UHFFFAOYSA-N").get("Synonyms").get(0));
		
		
//		ArrayList<Biotransformation> bt = hsbt.simulateOneStep(molecule);
//		System.out.println("Number of biotransformations: " + bt.size());
//		hsbt.cyb.saveBioTransformationProductsToSdf(bt, "data/BTransformerPaperTest/PA(10_0-a-13_0)_simulatorTest.sdf");
//		hsbt.simulateMetabolismAndSaveToSDF(molecule, "data/INRIA_Set/JarleiSet.sdf");
//		hsbt.simulateHumanMetabolismFromSDF("data/alpha-hydroquinone.mol","data/");
		
//		IAtomContainerSet metabolites = FileUtils.parseSdf("data/Curcumin_simulatorTest.sdf");
//		BufferedWriter bw0 = new BufferedWriter(new FileWriter("data/Bortezomib_simulatorTest.tsv"));		
//		for(IAtomContainer a : metabolites.atomContainers()){
//			bw0.write(a.getProperty("InChI") + "\t" + a.getProperty("InChIKey") + "\t" +
//					a.getProperty("Major Isotope Mass") + "\t" + a.getProperty("\t") + "\t" + 
//					a.getProperty("Metabolite ID") + "\t" + a.getProperty("\t") + "\t" + 
//					a.getProperty("Precursor InChIKey"));
//			bw0.write("\n");
//		}
//		bw0.close();
		
		
//		for(MetabolicReaction r : hsbt.cyb.reactionsByGroups.get("cypReactions")){
//			System.out.println(r.name);
//		}
//		System.out.println("N_DEALKYLATION_OF_PHOSPHORAMIDE: " +  ChemStructureExplorer.compoundMatchesReactionConstraints(hsbt.cyb.reactionsHash.get(ReactionName.N_DEALKYLATION_OF_PHOSPHORAMIDE.toString()),mt));
//		IAtomContainerSet acMetabolites = hsbt.cyb.generateAllMetabolitesFromAtomContainer(mt, 
//				hsbt.cyb.reactionsHash.get(ReactionName.N_DEALKYLATION_OF_PHOSPHORAMIDE.toString()), true);

		
//		BufferedReader bReader = new BufferedReader(new FileReader("/Users/yandj/Programming/Projects/SpectraPrediction/data/CASMI2016_original/CASMI2016-Cat3-Test_edit.tsv"));
//		BufferedWriter bWriter = new BufferedWriter(new FileWriter("/Users/yandj/Programming/Projects/SpectraPrediction/data/CASMI2016_original/CASMI2016-Cat3-Test_edit_with_mass.tsv"));
//		
//		String line;
//		
//		while((line = bReader.readLine()) != null ){
//			
//			String[] sline = line.split("\t");
//			IAtomContainer molecule = hsbt.cyb.getSmiParser().parseSmiles(sline[1]);
//			Double mass = ChemStructureExplorer.getMajorIsotopeMass(molecule);
//			
//			bWriter.write(sline[0] + "\t" + sline[1] + "\t" + mass.toString());
//			bWriter.newLine();
//		}
//		
//		bWriter.close();
//		bReader.close();
		
	}

}
