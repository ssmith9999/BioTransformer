package biotransformer.biomolecule;

import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ambit2.smarts.SMIRKSManager;
import ambit2.smarts.SMIRKSReaction;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.utils.ChemStructureManipulator;

public class smirksTest {

	public smirksTest() {
		// TODO Auto-generated constructor stub	
	}
	
	public static void main(String[]args) throws Exception{
		IChemObjectBuilder 	builder 	= SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 			= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		SMIRKSManager smrkMan = new SMIRKSManager(SilentChemObjectBuilder.getInstance());
		MDLV2000Writer mw = new MDLV2000Writer(System.out);
		
		IAtomContainer mol = smiParser.parseSmiles("OC[C@@H]1C[C@H](O)[C@@H](O)C(OC2=CC(OS([O-])(=O)=O)=C([N-][N+]#N)C(OP([O-])([O-])=O)=C2)O1");
//		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
		SMIRKSReaction sr = smrkMan.parse("[#8;X1-:6][P;X4:2]([#8;X1-:1])(=[O;X1:3])[#8;X2:4]-[*,#1:5]>>[H][#8;X2:1][P;X4:2](=[O;X1:3])([#8;X2:6][H])[#8;X2:4]-[*,#1:5]");
		IAtomContainerSet transforms = smrkMan.applyTransformationWithCombinedOverlappedPos(mol, null, sr);

		int i = 0;
		for(IAtomContainer atc : transforms.atomContainers()){
			i++;
			System.out.println(i);
//			System.out.println(smiGen.create(atc));
			
//			The SMILES gneration provides the following error, no matter whether I perceive the atom types or not.
//			Exception in thread "main" java.lang.NullPointerException: One or more atoms had an undefined number of implicit hydrogens
						
			mw.writeMolecule(atc);
		}
		
		// New adding explicit hydrogens
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
		IAtomContainerSet transforms2 = smrkMan.applyTransformationWithCombinedOverlappedPos(mol, null, sr);

		for(IAtomContainer atc : transforms.atomContainers()){
			System.out.println(i);
//			System.out.println(smiGen.create(atc));
			mw.writeMolecule(atc);
		}
		
		IAtomContainer mol2 = smiParser.parseSmiles("OC[C@@H]1C[C@H](O)[C@@H](O)C(OC2=CC(OS([O-])(=O)=O)=C([N-][N+]#N)C(OP([O-])([O-])=O)=C2)O1");
		Biotransformer bt = new Biotransformer(BioSystemName.HUMAN);
		IAtomContainer stmol = ChemStructureManipulator.standardizeMoleculeWithCopy(mol2);
		System.out.println(smiGen.create(stmol));
	}

}
