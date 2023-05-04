package in.sunfox.healthcare.java.core;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;

public class ParameterCalculations {
    public static void main(String[] args) throws IOException {

        String path = "/Users/sunfoxtechnologiessw/desktop/Files/";
        for (int i = 140; i <= 178; i++) {
            ArrayList<Double> lead2 = signalRead(path + "" + i + "/lead2_data.txt");
            EcgCharacteristicsCalculator ecgCharacteristicsCalculator=new EcgCharacteristicsCalculator(lead2,7,true,false);



//                EcgProcess ecgProcessor=new EcgProcess(ProcessorType.LEAD_TWO, new LeadTwoData(lead2),true,false);

//                EcgProcessorResult r= ecgProcessor.process(null,null,25);
            var ecgCharacteristics=ecgCharacteristicsCalculator.getCharacteristics(10000.0 / (lead2.size()-8));

            assert ecgCharacteristics != null;
            ArrayList<Object> result  =AshishWalia.characteristics(ecgCharacteristics);
            System.out.println(result);


        }
    }


    public static ArrayList<Double> signalRead(String signalPath) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(signalPath));
        ArrayList<Double> signalReadList = new ArrayList<Double>();

        String str = br.readLine();
        String[] arrayOfElements = str.split(",");

        if (arrayOfElements.length > 1) {
            for (String arrayOfElement : arrayOfElements) {
                signalReadList.add(Double.parseDouble(arrayOfElement));
            }
        } else {
            do {
                signalReadList.add(Double.parseDouble(str));
            } while ((str = br.readLine()) != null);
        }
        return signalReadList;

    }



}
