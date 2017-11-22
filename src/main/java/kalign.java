import model.SequenceInfo;
import util.FileHelper;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
/*Comment 1*/

/**
 * Created by nammi on 22/11/17.
 */
public class kalign {

    public static int numseq;
    public static int numprofiles;

    public static void main(String[] args) {
        SequenceInfo si = new SequenceInfo();
        si = readSequence("input.fasta", si);

    }

    public static SequenceInfo readSequence(String queryFile, SequenceInfo si){

        HashMap<String,String> sequences = new HashMap<String, String>();
        sequences = FileHelper.readSequenceFile(queryFile, sequences);
        Iterator iterator =  sequences.entrySet().iterator();
        int i=0;
        numseq = sequences.size();
        numprofiles = (numseq *2 ) - 1;
        si.s =  new String[sequences.size()];
        si.sn =  new String[sequences.size()];
        si.sl = new int[numprofiles];
        si.lsn =  new int[numseq];
        while (iterator.hasNext()){
            Map.Entry pair = (Map.Entry) iterator.next();
            //TODO : Extra last space not allocated
            String key = (String) pair.getKey();
            String value = (String) pair.getValue();
            si.s[i] = key;
            si.sl[i] = key.length();
            si.sn[i] = value;
            si.lsn[i] = value.length();
            i++;
        }

        int c = 0;
        int n = 0;
        i = 0;
        int j = 0;



        si.sip = new int[numprofiles][];
        si.nsip = new int[numprofiles];
        si.gis = new int[numprofiles][];
        si.relpos = new int[numprofiles][];

        for(i= numseq-1 ;i>=0; i--) {
            //TODO : Extra last space allocated
            si.gis[i] = new int[si.sl[i] + 1];
            si.relpos[i] = new int[si.sl[i] + 1];


            for (j = si.sl[i]; j >= 0; j--) {
                si.gis[i][j] = 0;
                si.relpos[i][j] = j;
            }
        }
        return si;
    }
}


