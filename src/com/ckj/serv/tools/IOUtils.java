package com.ckj.serv.tools;

import java.io.*;

/**
 *
 * @author Andy.Chen
 * @mail Chenjunjun.ZJ@gmail.com
 *
 */
public class IOUtils {

    public static String ReadString(InputStream inputStream){
        try {
            BufferedReader reader=new BufferedReader(new InputStreamReader(inputStream,"utf-8"));
            StringBuilder builder=new StringBuilder();
            String line=null;
            while ((line=reader.readLine())!=null){
                builder.append(line);
            }
            return builder.toString();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }
}