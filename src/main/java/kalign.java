import model.SequenceInfo;
import util.FileHelper;
import util.MatrixHelper;

import java.util.*;
/*Comment 1*/

/**
 * Created by nammi on 22/11/17.
 */
public class kalign {

    public static int numseq;
    public static int numprofiles;
    public static int[][][] matches = new int[8000][][];

    public static void main(String[] args) {
        SequenceInfo si = new SequenceInfo();
        si = readSequence("in.fasta", si);
        for (int i = 8000-1;i>=0;i--){
            matches[i] = null;
        }
        fill_hash(si,matches);
    }

    public static SequenceInfo readSequence(String queryFile, SequenceInfo si){

        LinkedHashMap<String,String> sequences = new LinkedHashMap<String, String>();
        sequences = FileHelper.readSequenceFile(queryFile, sequences);
        Iterator iterator =  sequences.entrySet().iterator();
        int i=0;
        numseq = sequences.size();
        numprofiles = (numseq *2 ) - 1;
        si.s =  new int[numseq][];
        si.sn =  new String[sequences.size()];
        si.sl = new int[numprofiles];
        si.lsn =  new int[numseq];
        while (iterator.hasNext()){
            Map.Entry pair = (Map.Entry) iterator.next();
            //TODO : Extra last space not allocated
            String key = (String) pair.getKey();
            String value = (String) pair.getValue();
            String seq =  key.toUpperCase();
            // TODO : Allocating 1 more letter space
            si.s[i] = new int[key.length()+1];
            for (int j=0;j<key.length();j++){
                si.s[i][j] = seq.charAt(j)-65;
            }
            // TODO : assigning 0 to last extra space
            si.s[i][key.length()] = 0;
            MatrixHelper.printArray(si.s[i]);
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

    public static void fill_hash(SequenceInfo si,int matches[][][])
    {
        int i,j,c,f;
        int key = 0;
        // size 26
	    int aacode[] = {0,0,1,2,3,4,5,6,7,-1,8,9,10,11,-1,12,13,14,15,16,-1,17,18,0,19,0};
        // size 20
        int aadecode[] = {0,2,3,4,5,6,7,8,10,11,12,13,15,16,17,18,19,21,22,24};
        //size 30
        int patterns[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        // size 10
	    int map[] = {0,0,0,0,0,0,0,0,0,0};
        int[] sp;
        int out[][] = null;
        System.out.println("Printing Keys\n");
        for (i = numseq-1;i>=0;i--){
            sp = si.s[i];
            for(j = si.sl[i]-2-1;j>=0;j--){
                key = ((aacode[sp[j]]*20)+aacode[sp[j+1]])*20+aacode[sp[j+2]];
                //System.out.print(key+" ");
                if (matches[key]==null){
                   // System.out.println("key "+ key);
                    matches[key] = new int [numprofiles][];
                    for (c = numprofiles-1;c>=0;c--){
                        matches[key][c] = null;
                    }
                }
                if (matches[key][i] == null){
                    System.out.println("key " + key + " i " + i + " j " + j);
                    matches[key][i] = new int[10];
                    matches[key][i][0] = 1;
                }
                if ((matches[key][i][0] %10) == 0){
                   System.out.println("Realloc key " + key + " i " + i + " j " + j);
                    matches[key][i] = Arrays.copyOf(matches[key][i], (matches[key][i][0] + 10));
                }
                matches[key][i][matches[key][i][0]] = j;
                matches[key][i][0] += 1;
            }
            System.out.println("\n");
        }

        for (i = 8000-1;i>=0;i--){
            if (matches[i] !=null){
                for (j = numseq-1;j>=0;j--){
                    if (matches[i][j] != null){
                        matches[i][j][0] |= 0x00040000;
                    }
                }
            }
        }

        c = 0;
        for (i = numseq-1;i>=0;i--){
            for (j = 8000-1;j>=0;j--){
                if (matches[j] != null){
                   // System.out.println("out j " + j);
                    if(matches[j][i] == null){
                       // System.out.println("j " + j + " i " + i +"\n");
                        map[c] = j;
                        key = j;
                        patterns[c*3 + 2] = aadecode[key / 400];
                        key %= 400;
                        patterns[c*3 + 1] = aadecode[key /20];
                        key %= 20;
                        patterns[c*3 + 0] = aadecode[key];
                        c++;
                    }
                    if(c == 10){
                        //start of 10x Wu-Manber;
                        out = ten_wu_manber(si.s[i],si.sl[i],patterns);
                        System.out.println("c = 10 ");
                        MatrixHelper.printMatrix(out);
                        for (f = 0;f < 10;f++){
                            matches[map[f]][i] = out[f];
                            matches[map[f]][i][0] |= 0x00010000;
                        }
                        //free(out);
                        c = 0;
                    }
                }
            }
            if (c!=0){
                for (f = c*3;f < 30;f++){
                    patterns[f] = 9;
                }
                out = ten_wu_manber(si.s[i],si.sl[i],patterns);
                System.out.println("c!=0");
                MatrixHelper.printMatrix(out);
                for (f = 0;f < c;f++){
                    matches[map[f]][i] = out[f];
                    matches[map[f]][i][0] |= 0x00010000;
                }
                //free(out);
                c = 0;
            }
        }
        MatrixHelper.printArray(patterns);

    }

    public static int[][] ten_wu_manber(int[] seq,int len,int p[])
    {
        // unsigned
        int Tc = 0;
        int Tc2 = 0;
        // register
        int R0 = 0;
        int R1= 153391689;//1;
        int S0 = 0;
        int S1 = 0;
        // size 26
        int T[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        // size 26
        int T2[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	    int [][]out = null;
        out = new int[10][];
        for (Tc = 10-1;Tc>=0;Tc--){
            out[Tc] = new int[10];
            out[Tc][0] = 1;
        }

        for (Tc = 0; Tc < 30;Tc++){
            T[p[Tc]] |= 1 << Tc;
            T2[p[Tc]] = T[p[Tc]]&153391689;//MASK;
        }
        while(len>0){
            //fprintf(stderr,"%d\n",len);
            len--;
            Tc = T[seq[len]];
            Tc2 = T2[seq[len]];
            S0 <<= 1;
            S0 |= 153391689;//MASK;
            S0 &= Tc;
            S0 |= Tc2;
            //S1 = ((R1 << 1) & Tc) | ((R0 | S0) << 1) | R0;//| 1;
            S1 = (R1 << 1) & Tc;//| 1;
            S1 |= (R0 | S0) << 1;
            S1 |= R0;
            S1 |= 153391689;

            //S1 = (R1 << 1) & (Tc) | ((R0 | S0) << 1) | (R0) | (153391689);

            if ((S1 & 0x24924924) !=0 ){
                if ((S1 & 0x00004924) != 0 ){
                    if ((S1 & 0x00000024) !=0){
                        if ((S1 & 4) !=0){
                            if ((out[0][0] %10) !=0){
                                out[0] = Arrays.copyOf(out[0],(out[0][0]+10));
                            }
                            out[0][out[0][0]] = len;
                            out[0][0]++;
                        }
                        if ((S1 & 32) !=0){
                            if ((out[1][0] %10) !=0){
                                out[1] = Arrays.copyOf(out[1],(out[1][0]+10));
                            }
                            out[1][out[1][0]] = len;
                            out[1][0]++;
                        }
                    }
                    if ((S1 & 0x00004900) !=0){
                        if ((S1 & 256)!=0){
                            if ((out[2][0] %10) !=0){
                                out[2] = Arrays.copyOf(out[2],(out[2][0]+10));
                            }
                            out[2][out[2][0]] = len;
                            out[2][0]++;
                        }
                        if ((S1 & 2048) !=0){
                            if ((out[3][0] %10) !=0){
                                out[3] = Arrays.copyOf(out[3],(out[3][0]+10));
                            }
                            out[3][out[3][0]] = len;
                            out[3][0]++;
                        }
                        if ((S1 & 16384)!=0){
                            if ((out[4][0] %10) !=0){
                                out[4] = Arrays.copyOf(out[4],(out[4][0]+10));
                            }
                            out[4][out[4][0]] = len;
                            out[4][0]++;
                        }
                    }
                }



                if ((S1 & 0x24920000) !=0){
                    if ((S1 & 0x00920000) !=0){
                        if ((S1 & 131072) !=0){
                            if ((out[5][0] %10) !=0){
                                out[5] = Arrays.copyOf(out[5],(out[5][0]+10));
                            }
                            out[5][out[5][0]] = len;
                            out[5][0]++;
                        }
                        if ((S1 & 1048576) !=0){
                            if ((out[6][0] %10) !=0){
                                out[6] = Arrays.copyOf(out[6],(out[6][0]+10));
                            }
                            out[6][out[6][0]] = len;
                            out[6][0]++;
                        }
                        if ((S1 & 8388608) !=0){
                            if ((out[7][0] %10) !=0){
                                out[7] = Arrays.copyOf(out[7],(out[7][0]+10));
                            }
                            out[7][out[7][0]] = len;
                            out[7][0]++;
                        }
                    }

                    if ((S1 & 0x24000000) !=0){
                        if ((S1 & 67108864) !=0){
                            if ((out[8][0] %10) !=0){
                                out[8] = Arrays.copyOf(out[8],(out[8][0]+10));
                            }
                            out[8][out[8][0]] = len;
                            out[8][0]++;
                        }
                        if ((S1 & 536870912) !=0){
                            if ((out[9][0] %10) !=0){
                                out[9] = Arrays.copyOf(out[9],(out[9][0]+10));
                            }
                            out[9][out[9][0]] = len;
                            out[9][0]++;
                        }
                    }
                }
            }
            R0 = S0;
            R1 = S1;
        }
        return out;
    }
}


