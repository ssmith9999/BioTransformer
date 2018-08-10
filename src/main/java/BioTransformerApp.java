import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.MissingOptionException;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;
import biotransformer.transformation.Biotransformation;
import biotransformer.utils.FileUtilities;
import biotransformer.utils.HumanSuperBioTransformer;
import biotransformer.utils.MetaboliteFinder;
import biotransformer.utils.MetaboliteFinder.FinderOption;

import org.apache.commons.io.FilenameUtils;


public class BioTransformerApp {

	
	public BioTransformerApp() {
		// TODO Auto-generated constructor stub
		
	}
	
	private static Options generateOptions(){
		
		final Option taskOption = Option.builder("k")
				.required(true)
				.hasArg(true)
				.argName("Task")
				.longOpt("task")
				.desc("The task to be permed: pred for prediction, or cid for compound identification ")
				.build();
		
		
		final Option biotransformerOption = Option.builder("b")
				.required(true)
				.hasArg(true)
				.argName("BioTransformer Type")
				.longOpt("btType")
				.desc("The type of description: Type of biotransformer - EC-based  (ecbased), CYP450 (cyp450), Phase II (phaseII), "
						+ "Human gut microbial (hgut), human super transformer* (superbio, or allHuman), Environmental microbial (envimicro)**.\n"
						+ "If option -m is enabled, the ovly valid biotransformer types are allHuman and env.")
				.build();

		final Option formatOption = Option.builder("f")
				.required(true)
				.hasArg(true)
				.argName("Input format")
				.longOpt("format")
				.desc("The format of the input: SMILES (smi), MDL Mol file (mol), or Structure data file (sdf).")
				.build();
		
		final Option nrOfStepsOption = Option.builder("s")
				.required(false)
				.hasArg(true)
				.argName("Number of steps")
				.longOpt("nsteps")
				.desc("The number of steps for the prediction. This option can be set by the user for the EC-based, CYP450, Phase II, and Environmental microbial biotransformers. The default value is 1.")
				.build();
		
		final Option inputOption = Option.builder("i")
				.required(true)
				.hasArg(true)
				.argName("Input")
				.longOpt("input")
				.desc("The input, which can be a SMILES string, a Mol file, or SDF file.")
				.build();
		
		final Option outputOption = Option.builder("o")
				.required(false)
				.hasArg(true)
				.argName("Output")
				.longOpt("ioutput")
				.desc("The output file name (which must be a SDF file) for single input queries, or the output folder for multiple input queries.\n"
						+ "When submitting a mutiple input query, an output folder must be provided.")
				.build();
		

		final Option annotateOption = Option.builder("a")
				.required(false)
				.hasArg(false)
				.argName("annotate")
				.longOpt("annotate")
				.desc("Search PuChem for each product, and annotate with CID and synonyms.")
				.build();
		
		final Option indentificationMassOption = Option.builder("m")
				.required(false)
				.hasArg(true)
				.argName("masses")
				.longOpt("masses")
				.desc("Semicolon-separated list of masses of compounds to identify")
				.build();
		
		final Option indentificationFormulaMetadataOption = Option.builder("r")
				.required(false)
				.hasArg(true)
				.argName("formulas")
				.longOpt("formulas")
				.desc("Semicolon-separated list of formulas of compounds to identify")
				.build();
		
//		final Option indentificationMetadataOption = Option.builder("d")
//				.required(false)
//				.hasArg(true)
//				.argName("Metadata used for identification")
//				.longOpt("criteria")
//				.desc("The metadata that is used for identification. It must be wither 'mass', or 'formula', or 'combo' ")
//				.build();		
//		
//		final Option formulaOption = Option.builder("r")
//				.required(false)
//				.hasArg(true)
//				.argName("A semicolon-separated list of masses and/or formulas to identify")
//				.longOpt("dList")
//				.desc("Given the starting comound(s), Find all their metabolites with the specified masses, formulas, or combinations thereof, and show a metabolic pathway. "
//						+ "A semicolon-separated list of masses is expected. E.g. '308.08:C15H16O7;', 'C6H6O2;C15H16O7', '308.08:320.09'"
//						+ "If a combination of mass and formula is given, they shoulld be separated by a semicolon.")
//				.build();	
//		
//		final Option massformulaOption = Option.builder("l")
//				.required(false)
//				.hasArg(true)
//				.argName("List of masse/formulas combinations to identify")
//				.longOpt("mList")
//				.desc("Given the starting comound(s), Find all their metabolites with the specified masses and formulas, and show a metabolic pathway. A semicolon-separated list of masses is expected. E.g. '308.08:C15H16O7;C6H6O2;320.09'")
//				.build();			
		
		final Option massToleranceOption = Option.builder("t")
				.required(false)
				.hasArg(true)
				.argName("Mass Tolerance")
				.longOpt("mTolerance")
				.desc("Mass tolerance for metabolite identification (default is 0.0).")
				.build();

		final Option helpOption = Option.builder("h")
				.required(false)
				.hasArg(false)
				.argName("help")
				.longOpt("help")
				.desc("Prints the usage.")
				.build();
		
		final Options options = new Options();
		options.addOption(taskOption);
		options.addOption(biotransformerOption);
		options.addOption(formatOption);
		options.addOption(nrOfStepsOption);
		options.addOption(inputOption);
		options.addOption(outputOption);
		options.addOption(indentificationMassOption);
		options.addOption(indentificationFormulaMetadataOption);	
//		options.addOption(massOption);
//		options.addOption(formulaOption);
//		options.addOption(massformulaOption);
		options.addOption(massToleranceOption);
		options.addOption(annotateOption);
		options.addOption(helpOption);

		return options;
	}
	
	public static CommandLine generateCommandLine(
			final Options options, final String[] commandLineArguments) throws ParseException{
		final CommandLineParser cmdLineParser = new DefaultParser();
		CommandLine commandLine = null;
		
		String header = "\nThis is the version 1.0.5 of BioTransformer. BioTransformer is a software tool that predicts small molecule metabolism in mammals, their gut microbiota,"
				+ " as well as the soil/aquatic microbiota. BioTransformer also assists scientists in metabolite identification, based on the metabolism prediction. \n\n";
		
		String footer = "* While the superbio option runs a set number of transformation steps in a pre-defined order (e.g. deconjugation first, then Oxidation/reduction, etc.), "
				+ "the allHuman option predicts all possible metabolites from any applicable reaction(Oxidation, reduction, (de-)conjudation) at each step."
				+ "\n\n** For the environmental microbial biodegradation, all reactions (aerobic and anaerobic) are reported, and not only the aerobic biotransformations "
				+ "(as per default in the EAWAG BBD/PPS system)."
				+ "\n\n*********\nExamples:\n*********\n\n1) To predict the biotransformation of a molecule from an SDF input using the human super transformer, use java -jar biotransformer-1-0-5.jar -k pred -b superbio "
				+ "-f sdf -i #{input file name} -o #{output folder}.\n\n"
				+ "\t\t- For each of the query molecule in the input file, an outputfile will be created with the list of corresponding metabolites.\n\n"
				+ "2) To predict the 2-step biotransformation of Thymol (a monoterpene) using the human super transformer (option allHuman) using the SMILES input, run\n\n"
				+ "java -jar biotransformer-1-0-5.jar -k pred -b allHuman -f smiles -i \"CC(C)C1=CC=C(C)C=C1O\" -o #{replace with output file name} -s 2\n\n"
				+ "Currently, the outputfile is SDF per default.\n\n"
				+ "3) Identify all human metabolites (max depth = 2) of Epicatechin (\"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\") with masses 292.0946 Da and 304.0946 Da, with a mass tolerance "
				+ "of 0.01 Da. Provide an annotation (Common name).\n\n"
				+ "java -jar biotransformer-1-0-5.jar -k cid -b allHuman -f smiles -i \"O[C@@H]1CC2=C(O)C=C(O)C=C2O[C@@H]1C1=CC=C(O)C(O)=C1\" -o #{replace with output file name} "
				+ "-s 2 -m \"292.0946; 304.0946\" -t 0.01 -a\n\t\t- DO NOT forget the quotes around the SMILES string or the list of masses.\n\n"
				+ "To report issues, provide feedback, or ask questions, please send an e-mail the following address: djoumbou@ualberta.ca\n\n"
				+ "BioTransformer is offered to the public as a freely acessible software package. Beside the prediction software, a manually curated database called BioTransformerDB is also available.\n\n"
				+ "Users are free to copy and redistribute the material in any medium or format. Moreover, they could modify, and build upon the material unfer the condition that they must give appropriate "
				+ "credit, provide links to the license, and indicate if changes were made. Furthermore, use and re-distribution of the these resources, in whole or in part, for commercial purposes requires "
				+ "explicit permission of the authors. We ask that all users of the BioTransformer software tool or BioTransformerDB to cite the BioTransformer reference in any resulting publications, and to acknowledge the authors.\n\n"
				+ "Cite: Djoumbou Feunang, Yannick; Cheminformatics Tools for Enabling Metabolomics; 2017; PhD Thesis";
		
		HelpFormatter formatter = new HelpFormatter();

		
		try{
			commandLine = cmdLineParser.parse(options, commandLineArguments);
		}
		catch (MissingOptionException missingOptionException){
			
			if( Arrays.asList(commandLineArguments).contains("-h") || Arrays.asList(commandLineArguments).contains("--help")){
				formatter.printHelp("\njava -jar biotransformer-1.0.5", header, options, footer, true);
			}
			else {
				System.out.println(missingOptionException.getLocalizedMessage());
			}			
		}
		catch (ParseException parseException){
			System.out.println("Could not parse the command line arguments "
					+ Arrays.toString(commandLineArguments) + "\nfor the following reaons:" 
					+ parseException);		
//			throw (parseException);
		}
		return commandLine;
	}
	
	
	
	public static void main(String[] args) throws Exception{
		IChemObjectBuilder 	builder = SilentChemObjectBuilder.getInstance();
		SmilesParser	smiParser		= new SmilesParser(builder);

		Options options = generateOptions();
		CommandLine commandLine = generateCommandLine(options, args);

//		System.out.println(commandLine.getOptionValue("b"));
//		System.out.println(commandLine.getOptionValue("f").length());
//		System.out.println(commandLine.getOptionValue("m"));
		
		IAtomContainer singleInput = null;
		String inputFileName = null;
//		IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		int nrOfSteps = 1;
		
		boolean annotate = false;
		String masses = null;
		String formulas = null;
		String metadata_input = null;		
		Double massTolerance = 0.01;
		String task = null;
		FinderOption opt = null;
		
		if(Arrays.asList(args).contains("-a") || Arrays.asList(args)
				.contains("--annotate")){		
			annotate = true;
			
			System.out.println("\n\n=============================================================");
			System.out.println("Compounds will be annotated using the PubChem API. Make sure");
			System.out.println("to have a secure internet connection.");
			System.out.println("=============================================================\n");
		}
		
		if(commandLine !=null){
			
			
			if(commandLine.getOptionValue("k") != null){
				task = commandLine.getOptionValue("k").trim();
				if( !(task.contentEquals("pred") || task.contentEquals("cid")) ){
					throw new IllegalArgumentException("You entered an invalid task. Enter either 'pred' (prediction) or 'cid' (compound identification)");
				}				
			}
			else {
				throw new IllegalArgumentException("Task type is missing. You must select either 'pred' (prediction) or 'cid' (compound identification)");
			}
				
//			System.out.println("OPTIONS");
			if(commandLine.getOptionValue("m") != null){
				masses = commandLine.getOptionValue("m").trim();
				if(masses.length() == 0){
					throw new IllegalArgumentException("You did not enter any mass.");
				}			
//				System.out.println("MASS: " + masses);
			}
			
			if(commandLine.getOptionValue("r") != null){
				formulas = commandLine.getOptionValue("r").trim();
				if(formulas.length() == 0){
					throw new IllegalArgumentException("You did not enter any formula.");
				}			
//				System.out.println("MASS: " + masses);
			}
			

			if(task == "cid"){
				if(masses != null && formulas != null){
					throw new IllegalArgumentException("You must enter either masses or formulas, not both");
				}	
				else if(masses == null && formulas == null){
					throw new IllegalArgumentException("Identification metadata are missing. Please add a list of masses (-m) or a list of formulas (-r)");
				}
				else if(masses != null){
					opt = FinderOption.MASS;
					metadata_input = masses;
				}
				else if(formulas != null){
					opt = FinderOption.FORMULA;
					metadata_input = formulas;
				}
			}
			

			if(commandLine.getOptionValue("s") != null){
				nrOfSteps = Integer.valueOf(commandLine.getOptionValue("s"));
//				System.out.println("nrOfSteps: " + nrOfSteps);
			}
			

			
			if(commandLine.getOptionValue("t") != null){
				if(commandLine.getOptionValue("t").trim().length() == 0){
					throw new IllegalArgumentException("You did not enter any mass tolerance. The parameter will be set to the default value (0.0)");
				}
				else{
					massTolerance = Double.valueOf(commandLine.getOptionValue("t").trim());	
//					System.out.println("massTolerance: " + massTolerance);
				}
			}

			
			final String biotransformerType = commandLine.getOptionValue("b");		
			final String format = commandLine.getOptionValue("f").toString().trim();
			final String outputF = commandLine.getOptionValue("o");
			
			if(format.contentEquals("smiles")){
				String smi = commandLine.getOptionValue("i");
				singleInput = smiParser.parseSmiles(smi);			
			}
			else if(format.contentEquals("mol") || format.contentEquals("sdf")){
				inputFileName = commandLine.getOptionValue("i");
//				containers = FileUtils.parseSdf(inputFileName);
				if(outputF == null){
					throw new MissingOptionException("A destination folder must be specified when your query molecules are provided in a file (Molfile or SDF). For more information, type java -jar biotransformer-1.0.5 --help.");
				}
			}
			else {
				throw new IllegalArgumentException("You entered an invalid format. It must be one of 'smiles','mol', or 'sdf'. Type java -jar biotransformer-1.0.5 --help for help.");
			}
			
			
			if(task == "cid"){
				
				if(metadata_input != null){
//					System.out.println("IDENTIFICATION TASK");
					if(biotransformerType.contentEquals("allHuman")){
						
						String[] mArr = masses.trim().split(";");
						ArrayList<String> list_of_masses = new ArrayList<String>();
						
						for(int k = 0; k < mArr.length; k++){
							try{
//								System.out.println(Double.valueOf(mArr[k].trim()));
								list_of_masses.add(mArr[k].trim());
							}
							catch(Exception e){
								System.err.println(e.getMessage());
							}
						}
//						System.out.println("Number of masses: " + dmasses.size());
						
						if (singleInput !=null){
							MetaboliteFinder mtf = new MetaboliteFinder();
							mtf.findAllHumanMetabolites(singleInput, list_of_masses, massTolerance, nrOfSteps, annotate, outputF, opt);
						}
						else {
							MetaboliteFinder mtf = new MetaboliteFinder();
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
							
							for(IAtomContainer atc : containers.atomContainers()){
								containers.add(mtf.findAllHumanMetabolites(atc, list_of_masses, massTolerance, nrOfSteps, annotate, opt));
							}
							
							SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
							sdfWriter.write(metabolites);
							sdfWriter.close();
						}
						
					}
					else if(biotransformerType.contentEquals("env")){
						String[] mArr = masses.trim().split(" ");
						ArrayList<String> dmasses = new ArrayList<String>();
						
						for(int k = 0; k < mArr.length; k++){
							try{
								dmasses.add(mArr[k].trim());
							}
							catch(Exception e){
								System.err.println(e.getMessage());
							}
						}
						
						
						if (singleInput !=null){
							MetaboliteFinder mtf = new MetaboliteFinder();
							mtf.findAllEnvMicroMetabolites(singleInput, dmasses, massTolerance, nrOfSteps, annotate, outputF, opt);
						}
						else {
							MetaboliteFinder mtf = new MetaboliteFinder();
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
							
							for(IAtomContainer atc : containers.atomContainers()){
								containers.add(mtf.findAllEnvMicroMetabolites(atc, dmasses, massTolerance, nrOfSteps, annotate, opt));
							}
							
							SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
							sdfWriter.write(metabolites);
							sdfWriter.close();
						}					
					}

					else if(biotransformerType.contentEquals("superbio")){
						String[] mArr = masses.trim().split(" ");
						ArrayList<String> dmasses = new ArrayList<String>();
						
						for(int k = 0; k < mArr.length; k++){
							try{
								dmasses.add(mArr[k].trim());
							}
							catch(Exception e){
								System.err.println(e.getMessage());
							}
						}
							
						if (singleInput !=null){
							MetaboliteFinder mtf = new MetaboliteFinder();
							mtf.findSuperbioMetabolites(singleInput, dmasses, massTolerance, annotate, outputF, opt);
						}
						else {
							MetaboliteFinder mtf = new MetaboliteFinder();
							IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
							IAtomContainerSet metabolites = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
							
							for(IAtomContainer atc : containers.atomContainers()){
								containers.add(mtf.findSuperbioMetabolites(atc, dmasses, massTolerance, annotate, opt));
							}
							
							SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputF));		
							sdfWriter.write(metabolites);
							sdfWriter.close();
						}					
					}				
					else{
						throw new IllegalArgumentException("For metabolite identification, the biotransformer type must be either allHuman, superbio, or env.");
					}
					
								
					
				}
				else{
					throw new IllegalArgumentException("For metabolite identification, you must enter a list of masses, and/or molecular formulas");
				}
				
			}
			
			else if (biotransformerType.contentEquals("cyp450")){
				Cyp450BTransformer cyp450bt = new Cyp450BTransformer(BioSystemName.HUMAN);
				if (singleInput !=null){
					ArrayList<Biotransformation> biotransformations = cyp450bt.predictCyp450BiotransformationChain(singleInput, true, true, nrOfSteps, 0.5);
					
					if(outputF != null){
						cyp450bt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved in " + outputF + ".\n");
					}
					else{
						cyp450bt.saveBioTransformationProductsToSdf(biotransformations, "./query_CYP450_based_metabolites.sdf");
						System.out.println("Done!!!!\nThe result was saved to ./query_CYP450_based_metabolites.sdf");
					}
				}
				else {
					IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
						if (containers.getAtomContainerCount()>0){
						if (outputF	== null){
							throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
						}
						else{
							containers = FileUtilities.parseSdf(inputFileName);
							cyp450bt.simulateCyp450MetabolismAndSaveToSDF(containers, nrOfSteps, 0.5, outputF, annotate);
							System.out.println("Done!!!!\nThe result was saved into the specified folder.");
						}				
					}
				}
			}
			else if (biotransformerType.contentEquals("ecbased")){
				ECBasedBTransformer ecbt =  new ECBasedBTransformer(BioSystemName.HUMAN);
				
				if (singleInput !=null){				
					ArrayList<Biotransformation> biotransformations = ecbt.simulateECBasedMetabolismChain(singleInput, true, true, nrOfSteps, 0.5);				
					if(outputF != null){
						ecbt.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved in " + outputF + ".\n");
					}
					else{
						ecbt.saveBioTransformationProductsToSdf(biotransformations, "./query_EC_based_metabolites.sdf");
						System.out.println("Done!!!!\nThe result was saved to ./query_EC_based_metabolites.sdf");
					}				
				}
				else {
					IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
//					if (containers.getAtomContainerCount()>0) {
					System.err.println(containers.getAtomContainerCount() + " container(s).");
					if (outputF	!= null){
						ecbt.simulateECBasedMetabolismAndSaveToSDF(inputFileName, nrOfSteps, 0.5, outputF, annotate);	
						System.out.println("Done!!!!\nThe results were saved in " + outputF + ".\n");
					}
					else {
						throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
					}
				}
				
			}
			else if (biotransformerType.contentEquals("hgut")){
				HGutBTransformer hgut = new HGutBTransformer();
				if(nrOfSteps!=6){
					System.out.println("\n=======>The number of steps for reductive metabolism is fixed at 6. No need to set a number of steps for the human gut metabolism.\n\n");
				}
				
				if (singleInput !=null){				
					ArrayList<Biotransformation> biotransformations = hgut.simulateGutMicrobialMetabolism(singleInput, true, true, 8, 0.5);				
					if(outputF != null){
						hgut.saveBioTransformationProductsToSdf(biotransformations, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved in " + outputF + ".\n");
					}
					else{
						hgut.saveBioTransformationProductsToSdf(biotransformations, "./query_EC_based_metabolites.sdf", annotate);
						System.out.println("Done!!!!\nThe result was saved to ./query_EC_based_metabolites.sdf");
					}				
				}
				else {
//					if (containers.getAtomContainerCount()>0) {
//					System.err.println(containers.getAtomContainerCount() + " container(s).");
					if (outputF	!= null){
						
						hgut.simulateGutMicrobialMetabolismAndSave(inputFileName, true, 
								true, 8, 0.5, outputF + FilenameUtils.getBaseName(inputFileName).toString() + "_bioT_hgut_metabolites.sdf", annotate);
//						hgut.simulateGutMicrobialPolyphenolMetabolismAndSaveToSDF(containers, 6, 0.5, outputF);	
						System.out.println("Done!!!!\nThe results were saved in " + outputF + ".\n");
					}
					else {
						throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
					}
				}			
				
			}
			else if (biotransformerType.contentEquals("phaseII")){
				Phase2BTransformer phase2b = new Phase2BTransformer(BioSystemName.HUMAN);
				if (singleInput !=null){					
//					hsbt.simulateMetabolismAndSaveToSDF(singleInput, "./query_BioT_sim_metabolites.sdf");
					ArrayList<Biotransformation> biots =
					phase2b.applyPhase2TransformationsChainAndReturnBiotransformations(singleInput,
							true, true, true, nrOfSteps, 0.0);
					if(outputF != null){
						phase2b.saveBioTransformationProductsToSdf(biots, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved to " + outputF);

					}
					else
					{
						phase2b.saveBioTransformationProductsToSdf(biots, "./query_BioT_phaseII_metabolites.sdf", annotate);
						System.out.println("Done!!!!\nThe result was saved to ./query_BioT_phaseII_metabolites.sdf");

					}
				}
				else {
					IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
					
					if (containers.getAtomContainerCount()>0){
						if (outputF	== null){
							throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
						}
						else{
							phase2b.applyReactionChainFromSdfToSingleSdf(inputFileName, true, true, true, nrOfSteps, 0.0, outputF + FilenameUtils.getBaseName(inputFileName).toString() + "_phaseII_transformer_metabolites.sdf", annotate);
	//						hsbt.simulateHumanMetabolismAndSaveToSDF(containers, outputF);
							System.out.println("Done!!!!\nThe result was saved into the specified folder.");
						}				
					}
				}
			}
			else if (biotransformerType.contentEquals("superbio")){
				if(nrOfSteps!=6){
					System.out.println("\n\n=======>The configutration is fixed for this simulation. No need to set a number of steps for the super human transformer.\n\n");
				}

				HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();
				if (singleInput !=null){	
					
					if(outputF != null){
						hsbt.simulateHumanAndGutMicrobialMetabolismAndSaveToSDF(singleInput, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved to " + outputF);
						
					} else{
						hsbt.simulateHumanAndGutMicrobialMetabolismAndSaveToSDF(singleInput, "./query_BioT_sim_metabolites.sdf", annotate);
						System.out.println("Done!!!!\nThe result was saved to ./query_bioT_supertransformer_sim_metabolites.sdf");
						
					}
				}
				else {
//					IAtomContainerSet containers = FileUtils.parseSdf(inputFileName);
					
//					if (containers.getAtomContainerCount()>0){
						if (outputF	== null){
							throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
						}
						else{
							hsbt.simulateHumanMetabolismFromSDFtoSingleSDF(inputFileName, outputF + FilenameUtils.getBaseName(inputFileName).toString() + "_bioT_supertransformer_metabolites.sdf", annotate);
	//						hsbt.simulateHumanMetabolismAndSaveToSDF(containers, outputF);
							System.out.println("Done!!!!\nThe result was saved into the specified folder.");
						}				
					}
//				}
			}
			
			else if (biotransformerType.contentEquals("allHuman")){
				HumanSuperBioTransformer hsbt = new HumanSuperBioTransformer();

				if (singleInput !=null){
					
					if(outputF != null){
						hsbt.predictAllHumanBiotransformationChainAndSaveToSDF(singleInput, nrOfSteps, 0.5, outputF, annotate);
						System.out.println("Done!!!!\nThe result was saved in " + outputF + ".\n");
					}
					else{
						hsbt.predictAllHumanBiotransformationChainAndSaveToSDF(singleInput, nrOfSteps, 0.5, "./query_EC_based_metabolites.sdf", annotate);
						System.out.println("Done!!!!\nThe result was saved to ./query_EC_based_metabolites.sdf");
					}				
				}
				else {
//					if (containers.getAtomContainerCount()>0) {
//					System.err.println(containers.getAtomContainerCount() + " container(s).");
					if (outputF	!= null){
						hsbt.predictMetabolismAllHumanFromSDFAndSavetoSDF(inputFileName, outputF, nrOfSteps, 0.5, annotate);
						System.out.println("Done!!!!\nThe results were saved in " + outputF + ".\n");
					}
					else {
						throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
					}
				}

			}
			else if (biotransformerType.contentEquals("env")){
				EnvMicroBTransformer ebt = new EnvMicroBTransformer();
				if (singleInput !=null){
					ebt.simulateEnvMicrobialDegradationAndSaveToSDF(singleInput, true, true, nrOfSteps, 0.5, outputF, annotate);
					System.out.println("Done!!!!\nThe result was saved to " + outputF);
				}
				else {
					
					IAtomContainerSet containers = FileUtilities.parseSdf(inputFileName);
					
					if (containers.getAtomContainerCount()>0){
						if (outputF	== null){
							throw new MissingOptionException ("You must enter a destination folder when you submit a query in the following Molfile and SDF formats. Type java -jar biotransformer-1.0.5 --help for help.");
						}
						else{
							ebt.simulateEnvMicrobialDegradationAndSaveToSDF(containers, true, true, nrOfSteps, 0.5, outputF, annotate);
							System.out.println("Done!!!!\nThe result was saved into the specified folder.");
						}				
					}
			}
			}	
			else {
				throw new IllegalArgumentException("You entered an invalid biotransformer option.\n"
						+ "Choose one of the following: ecbased, cyp450, hgut, phaseII, superbio, allHuman, or env.\nType java -jar biotransformer-1.0.5 --help for help.");
			}			
		}
		
	}

}
