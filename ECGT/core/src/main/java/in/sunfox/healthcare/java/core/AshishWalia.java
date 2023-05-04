package in.sunfox.healthcare.java.core;

import java.util.ArrayList;

public class AshishWalia {
    public static ArrayList<Object> characteristics(EcgCharacteristics list)
    {
        ArrayList<Object> newlist=new ArrayList<Object>();
        newlist.add(list.getPr());
        newlist.add( list.getHeartRate());
        newlist.add(list.getPr());
        newlist.add(list.getQrs());
        newlist.add(list.getQt());
        newlist.add( list.getQtc());
        newlist.add(list.getStElevation());
        return newlist;
    }
}
