package biotransformr;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import biotransformer.fingerprint.ChemStructureFingerprinter;
import biotransformer.utils.MetaboliteFinder;

public class MetaboliteFinderTest extends MetaboliteFinder {

	
	public MetaboliteFinderTest() throws IOException, ParseException, CDKException {
		// TODO Auto-generated constructor stub
		super();
	}
	
	public static void main(String[] args) throws Exception{
		MetaboliteFinderTest mft = new MetaboliteFinderTest();
		IAtomContainer atc = mft.hsbt.smiParser.parseSmiles("O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1");
		ArrayList<String> masses = new ArrayList<String>();
		ArrayList<String> formulas = new ArrayList<String>();
//		masses.add("292.0946:C15H14O6");
//		masses.add("226.0841");
//		masses.add("166.0629");
//		masses.add("292.0946");
//		masses.add("192.0786");
//		formulas.add("C11H14O3");
//		masses.add("166.062375");
//		masses.add("241.1064");
//		masses.add("273.0758");
//		masses.add("315.0712");
//		masses.add("289.0708");
//		masses.add("289.0707");
//		masses.add("273.0427");
//		masses.add("305.1022");
//		masses.add("273.0428");
//		masses.add("193.0858");
//		masses.add("369.1181");
//		masses.add("177.0909");
//		masses.add("209.0808");
//		masses.add("403.1235");
//		masses.add("227.0914");
//		masses.add("385.1131");
//		masses.add("451.1234");
//		masses.add("209.0808");
//		masses.add("223.0963");
//		masses.add("305.1020");
//		masses.add("305.1019");
//		masses.add("196.0604");
//		masses.add("467.1184");

//		masses.add("110.0362");
//		masses.add("195.0526");
//		masses.add("179.0577");
//		masses.add("464.0638");
//		masses.add("326.0998");
//		masses.add("358.0949");
//		masses.add("208.0730");
//		masses.add("226.0836");
//		masses.add("402.1157");
//		masses.add("182.0572");
//		masses.add("222.0885");
//		masses.add("166.0624");
//		masses.add("342.0946");
//		masses.add("208.0730");
//		masses.add("384.1053");
//		masses.add("240.0993");
//		masses.add("260.1080");
//		masses.add("288.0629");
//		masses.add("368.1103");
//		masses.add("166.0624");
//		masses.add("192.0780");
//		masses.add("240.0986");
//		masses.add("210.0888");
//		masses.add("466.1106");
//		masses.add("272.0680");
//		masses.add("314.0634");
//		masses.add("122.0363");
//		masses.add("452.1282");
//		masses.add("288.0630");
//		masses.add("138.0313");
//		masses.add("450.1156");
//		masses.add("304.0942");
//		masses.add("272.0349");
//		masses.add("166.0625");
//		masses.add("178.0987");
//		masses.add("436.1346");
//		masses.add("304.0941");
//		masses.add("176.0831");
//		masses.add("224.1044");
//		masses.add("226.0837");
//		masses.add("208.0731");
//		masses.add("166.0626");
//		masses.add("272.0350");
//		masses.add("304.0944");
//		masses.add("258.1098");
//		masses.add("184.0731");
//		masses.add("480.1262");
//		masses.add("232.1309");
//		masses.add("464.0935");
		
		masses.add("110.0362");
		masses.add("122.0363");
		masses.add("138.0313");
		masses.add("166.0624");
		masses.add("176.0831");
		masses.add("178.0987");
		masses.add("179.0577");
		masses.add("182.0572");
		masses.add("184.0731");
		masses.add("192.0780");
		masses.add("195.0526");
		masses.add("208.0730");
		masses.add("210.0888");
		masses.add("222.0885");
		masses.add("224.1044");
		masses.add("226.0836");
		masses.add("232.1309");
		masses.add("240.0986");
		masses.add("258.1098");
		masses.add("260.1080");
		masses.add("272.0349");
		masses.add("272.0680");
		masses.add("288.0629");
		masses.add("304.0941");
		masses.add("314.0634");
		masses.add("326.0998");
		masses.add("342.0946");
		masses.add("358.0949");
		masses.add("368.1103");
		masses.add("384.1053");
		masses.add("402.1157");
		masses.add("436.1346");
		masses.add("450.1156");
		masses.add("452.1282");
		masses.add("464.0638");
		masses.add("464.0935");
		masses.add("466.1106");
		masses.add("480.1262");
		
		mft.findSuperbioMetabolites(atc, masses, 0.01, true, "/Users/yandj/epicatechin_metabolites_identification_38_001.sdf", FinderOption.MASS);
	
//		masses.add("208.0735");
//		masses.add("292.0946");
//		formulas.add("C11H12O3");
////		mft.findAllHumanMetabolites(atc, masses, 0.01, 2, false, "/Users/yandj/epicatechin_allHuman_met_id.sdf", FinderOption.MASS);
//		mft.findAllHumanMetabolites(atc, formulas, 0.01, 2, false, "/Users/yandj/epicatechin_allHuman_met_id.sdf", FinderOption.FORMULA);
////		mft.findSuperbioMetabolites(atc, masses, 0.01, false, "/Users/yandj/epicatechin_superbio_met_id.sdf", FinderOption.MASS);
////		mft.findSuperbioMetabolites(atc, masses, 0.01, true, "/Users/yandj/epicatechin_superbio_met_id.sdf", FinderOption.MASS);

//		masses.add("308.0896");
//		mft.findAllEnvMicroMetabolites(atc, masses, 0.01, 1, false, "/Users/yandj/epicatechin_env_micro_met_id.sdf", FinderOption.MASS);
//		formulas.add("C15H16O7");
//		mft.findAllEnvMicroMetabolites(atc, formulas, 0.01, 1, false, "/Users/yandj/epicatechin_env_micro_met_id.sdf", FinderOption.FORMULA);
//		ArrayList<String> masses_2 = new ArrayList<String>();
//		masses_2.add("308.08:C15H16O7");
////		mft.findAllEnvMicroMetabolites(atc, masses_2, 0.01, 1, false, "/Users/yandj/epicatechin_env_micro_met_id.sdf", FinderOption.MASSFORMULA);
////		ArrayList<String> formulas = new ArrayList<String>();
////		formulas.add("C15H14O6");
////		formulas.add();
////		mft.findAllHumanMetabolitesByMolecularFormula(atc, formulas, 2, false, "data/epicatechin_superbio_met_id_by_formula.sdf");
//		
//		IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
//		IMolecularFormula f1 = MolecularFormulaManipulator.getMolecularFormula("322.04", builder);
//		
//		System.out.println(f1 == null);
//		IMolecularFormula f2 = MolecularFormulaManipulator.getMolecularFormula("C6O2", builder);
//		System.out.println(f2 == null);
////		System.out.println(MolecularFormulaManipulator.getString(f1).contentEquals(MolecularFormulaManipulator.getString(f2)) );
//		System.out.println(MolecularFormulaManipulator.getString(f2));
////		System.out.println(MolecularFormulaManipulator.getString(f2));
//		
////		IAtomContainer atc2 = mft.hsbt.smiParser.parseSmiles("[H][C@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@]1([H])[C@]2([H])C=CC2=CC(=O)CC[C@@]12C");
////		LinkedHashMap<String,Integer> fp = ChemStructureFingerprinter.generateBTRailsCountFingerprint(atc2);
////		System.out.println(fp);
////		
////		IMolecularFormula mf = MolecularFormulaManipulator.getMolecularFormula("C6H2O", SilentChemObjectBuilder.getInstance());
////		System.out.println(mf);
////		System.out.println(MolecularFormulaManipulator.getString(mf));
		
	}

}
