/**
 * 
 * @author Djoumbou Feunang, Yannick, PhD
 *
 */

package biotransformer.utils;

import java.io.IOException;
import java.util.HashSet;
import java.util.LinkedHashMap;

import org.codehaus.jackson.JsonParseException;
import org.codehaus.jackson.map.JsonMappingException;
import org.json.simple.parser.ParseException;
import org.openscience.cdk.exception.CDKException;

import biotransformer.biomolecule.Enzyme;
import biotransformer.biosystems.BioSystem.BioSystemName;
import biotransformer.btransformers.Biotransformer;
import biotransformer.btransformers.Cyp450BTransformer;
import biotransformer.btransformers.ECBasedBTransformer;
import biotransformer.btransformers.EnvMicroBTransformer;
import biotransformer.btransformers.HGutBTransformer;
import biotransformer.btransformers.Phase2BTransformer;

public class Statistics {
	
	protected Biotransformer humanBT 		= new Biotransformer(BioSystemName.HUMAN);
	protected HGutBTransformer hgutBT 		= new HGutBTransformer();
	protected EnvMicroBTransformer envmBT 	= new EnvMicroBTransformer();
	protected Cyp450BTransformer hCyp450BT	= new Cyp450BTransformer(BioSystemName.HUMAN);
	protected Phase2BTransformer phaseII0BT	= new Phase2BTransformer(BioSystemName.HUMAN);
	protected ECBasedBTransformer ecBase450BT	= new ECBasedBTransformer(BioSystemName.HUMAN);

	public Statistics() throws JsonParseException, JsonMappingException, ParseException, IOException, CDKException{
		// TODO Auto-generated constructor stub
		generateStatistics();
	}
	
	public static LinkedHashMap<String, Integer> generalStatistics = new LinkedHashMap<String, Integer>();
	
	
	public void generateStatistics(){
		
		int totalNrOfEnzymesUsed = 0;
		int totalNrOfReactionsUsed = 0;
		int totalNumberOfPreferenceRules = 0;
		int totalNumberOfReactionsWithPreferenceRules = 0; // including the predominant and the predominated ones.
		
		
//		ArrayList<ReactionName> humanReactEC = 
		
//		generalStatistics.put("No. of enzymes mapped to the human biosystem ",humanBT.enzymesList.size());
		generalStatistics.put("No. of reactions mapped to the human biosystem ",humanBT.bSystem.getEnzymesList().size());
		generalStatistics.put("No. of enzymes mapped to the human gut microbiome biosystem ",hgutBT.bSystem.getEnzymesList().size());
		generalStatistics.put("No. of enzymes mapped to the environmental microbiome biosystem ",envmBT.bSystem.getEnzymesList().size());
//		System.out.println(hgutBT.bSystem.getEnzymesList());
		
		
		/*
		 * Metababolic pathways
		 */
		
		generalStatistics.put("No. of metabolic pathways mapped to the human biosystem", 
				humanBT.bSystem.getMetPathwaysHash().size());
		
//		int  metPathEnzAssociation = 0;
//		HashSet enzNames = new HashSet(); 
//		
////		for(Metabolic)
////		
////		generalStatistics.put("No. of metabolic pathways mapped to the human biosystem", 
////				humanBT.bSystem.getMetPathwaysHash().size());
////		metPathwaysHash

		System.out.println("\n\n"
		+ "/*******************\n"
		+ "      EC-based\n"
		+ "*******************/");
		ecBase450BT.printStatistics();

		
		System.out.println("\n\n"
				+ "/*******************\n"
				+ "      CYP450\n"
				+ "*******************/");
		hCyp450BT.printStatistics();
		
		System.out.println("\n\n"
				+ "/*******************\n"
				+ "      Human Gut\n"
				+ "*******************/");
		hgutBT.printStatistics();		
		
		System.out.println("\n\n"
				+ "/*******************\n"
				+ "      Phase II \n"
				+ "*******************/");
		phaseII0BT.printStatistics();
		
		System.out.println("\n\n"
				+ "/*******************\n"
				+ "     Environmental \n"
				+ "*******************/");
		envmBT.printStatistics();	
	}
	
	
	
	
//	public void

}
