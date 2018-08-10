package biotransformr;



import java.time.LocalDateTime;
import java.time.LocalTime;

import biotransformer.fingerprint.ChemStructureFingerprinter;

public class ChemStructureFingerPrinterTest extends ChemStructureFingerprinter{

	public ChemStructureFingerPrinterTest() {
		// TODO Auto-generated constructor stub
		super();
	}

	public static void main(String[] args) throws Exception{
		ChemStructureFingerPrinterTest csfp  = new ChemStructureFingerPrinterTest();
		long localTime_0 = System.currentTimeMillis();
		csfp.generateDEREP_NPFingerprint("/Users/yandj/Projects/DB_data/Drugbank/Drugbank-5.0/Drugbank-50.tsv");
		long localTime_1 = System.currentTimeMillis();
		
		System.out.println("Done after: " + (localTime_1 - localTime_0));
	
	}
	
}
