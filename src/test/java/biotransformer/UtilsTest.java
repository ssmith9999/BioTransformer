package biotransformer;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.utils.ChemStructureManipulator;;



public class UtilsTest {
	protected IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
	protected SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
	protected SmilesParser	smiParser		= new SmilesParser(builder);

	public UtilsTest(){
		
	}
	
	public static void main(String[] args) throws Exception{	
		UtilsTest utest = new UtilsTest();
		IAtomContainer molecule_1 = utest.smiParser.parseSmiles("CC([O-])=Nc1ccc(O)cc1");
		IAtomContainer molecule_2 = utest.smiParser.parseSmiles("CN(=O)=O");
		
		System.out.println(molecule_1);
		IAtomContainer smolecule_1 = ChemStructureManipulator.preprocessContainer(molecule_1);
		System.out.println(utest.smiGen.create(smolecule_1) );

		System.out.println(molecule_2);
		IAtomContainer smolecule_2 = ChemStructureManipulator.preprocessContainer(molecule_2);
		System.out.println(utest.smiGen.create(smolecule_2) );
		
	}
}
