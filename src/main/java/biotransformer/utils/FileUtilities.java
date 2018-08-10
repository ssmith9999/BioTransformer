package biotransformer.utils;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedHashMap;

import org.apache.commons.io.FilenameUtils;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;


public class FileUtilities {

	public FileUtilities() {
		// TODO Auto-generated constructor stub

	
	}
	
	
	public static IAtomContainerSet parseSdf(String sdfFileName) throws FileNotFoundException {
		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);
		
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(sdfFileName), bldr);
		
		while (sdfr.hasNext()){
			IAtomContainer mol = sdfr.next();
			containers.addAtomContainer(mol);	
		}
		
		
		return containers;
		
	}
	
	public static int countUniqueCompounds(String sdfFileName) throws FileNotFoundException, CDKException{
		int count = 0;
		
		IAtomContainerSet containers = parseSdf(sdfFileName);
		
		LinkedHashMap<String, IAtomContainer> lmh  = new LinkedHashMap<String, IAtomContainer>();
		
		for (IAtomContainer atc : containers.atomContainers()){
			String ikey = atc.getProperty("InChIKey");
			if(ikey == null){
				InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
				ikey = gen.getInchiKey();			
			}
			
			if(! lmh.containsKey(ikey)){
				lmh.put(ikey, atc);
				count++;
			}
		}
		
		
		return count;
	}

	
	public static void divideSdfFile(String sdfFileName, int limit) throws CDKException, IOException{
		IAtomContainerSet containers = parseSdf(sdfFileName);
		int part_nr = 1;
		int at_nr=0;
		
		if(containers.getAtomContainerCount()>limit){
			IAtomContainerSet molecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			System.out.println(sdfFileName);
			System.out.println(FilenameUtils.getBaseName(sdfFileName));
			System.out.println(FilenameUtils.getFullPathNoEndSeparator(sdfFileName));
			
			for(IAtomContainer atc : containers.atomContainers()){
				at_nr++;
				System.out.println(at_nr);
				molecules.addAtomContainer(atc);
				if(at_nr % limit == 0){
					SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(sdfFileName) + 
							"/" + FilenameUtils.getBaseName(sdfFileName)  + "_part_"+ part_nr + ".sdf" ));
					for(IAtomContainer a : molecules.atomContainers()){
						sdfWriter.write(a);						
					}
					sdfWriter.close();
					part_nr++;
					molecules.removeAllAtomContainers();
				}
			}
		
		
		}
	}


	public static void buildSdfFromTSV(String tsvFileName) throws IOException, CDKException{
		
		BufferedReader bRead = new BufferedReader(new FileReader(tsvFileName));
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(FilenameUtils.getFullPathNoEndSeparator(tsvFileName) + 
				"/" + FilenameUtils.getBaseName(tsvFileName)  + ".sdf" ));
		
		int counter = 0;
		String line = null;
		
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesGenerator smiGen 		= new SmilesGenerator().isomeric();
		SmilesParser	smiParser		= new SmilesParser(builder);
		InChIGeneratorFactory inchiGenFactory = InChIGeneratorFactory.getInstance();
		
//		System.out.println(counter);
		
		while((line = bRead.readLine()) !=null && !line.contains("SMILES")){
			counter++;
			System.out.println(counter);
			String[] sline = line.split("\t");
			
			if(!sline[1].contentEquals("NULL")){
				
				System.out.println(counter + " " + sline[0]);
				IAtomContainer atc= smiParser.parseSmiles(sline[1]);
				AtomContainerManipulator.suppressHydrogens(atc);
				InChIGenerator gen = inchiGenFactory.getInChIGenerator(atc);
				atc.setProperty(CDKConstants.TITLE, sline[0]);
				atc.setProperty("Name", sline[0]);
				
//				if(gen.getInchiKey().trim().contentEquals(sline[3].trim())){
					atc.setProperty("InChIKey", gen.getInchiKey());
//				}
//				else{
//					System.err.println("Issue of inchikey incompatibility with " + sline[0]);
//					System.err.println(gen.getInchiKey());
//					System.err.println(sline[0].trim());
//	//				break;
//				}
//				if(sline.length>=5){
//					atc.setProperty("Origin", sline[4]);
//				}
//				else{
//					atc.setProperty("Origin", null);
//				}
				
				atc.setProperty("SoMs", "");
				atc.setProperty("References","");
				sdfWriter.write( ChemStructureManipulator.preprocessContainer(atc) );
			}
		}
		
		sdfWriter.close();
	}

	
	public static void saveAtomContainersToSDF(IAtomContainerSet containers, String outputFileName) throws CDKException, IOException{
		SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));		
		sdfWriter.write(containers);
		sdfWriter.close();
	}

	public static String saveAtomContainersToString(IAtomContainerSet containers) throws CDKException, IOException{
		ByteArrayOutputStream containers_to_s = new ByteArrayOutputStream();
		SDFWriter sdfWriter = new SDFWriter(containers_to_s);		
		sdfWriter.write(containers);
		sdfWriter.close();
		
		return containers_to_s.toString();
		
	}
	

	public static String saveAtomContainerToString(IAtomContainer container) throws CDKException, IOException{
		ByteArrayOutputStream containers_to_s = new ByteArrayOutputStream();
		SDFWriter sdfWriter = new SDFWriter(containers_to_s);		
		sdfWriter.write(container);
		sdfWriter.close();
		
		return containers_to_s.toString();
		
	}	
}
