package com.ckj.serv.service;

import com.alibaba.fastjson.JSON;
import com.ckj.serv.model.ArrayReq;
import com.ckj.serv.model.ArrayResp;
import com.ckj.serv.tools.ArrayUitls;
import com.ckj.serv.tools.IOUtils;
import com.ckj.serv.tools.PcaUtils;
import org.ejml.data.DenseMatrix64F;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;

public class PcaService extends HttpServlet {

    @Override
    protected void service(HttpServletRequest req, HttpServletResponse resp) {
        try {
            String requestBody = IOUtils.ReadString(req.getInputStream());
            ArrayReq request = JSON.parseObject(requestBody, ArrayReq.class);
            double[][][] resp1 = pcaEncode(request);
            resp.getWriter().print(JSON.toJSONString(resp1));
            resp.getWriter().flush();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public double[][][] pcaEncode(ArrayReq arrayReq) {
        double[][][] result=new double[arrayReq.array.length][][];
        for(int i=0;i<arrayReq.array.length;i++){
            double[][] array2d = arrayReq.array[i];
            int columns = arrayReq.columns;
            DenseMatrix64F matrix64F = new DenseMatrix64F(array2d);
            DenseMatrix64F encArray = PcaUtils.pcaEncode(matrix64F, columns);
            double[] array = encArray.getData();
            double[][] enc2dArray = ArrayUitls.to2dArray(array, array2d.length);
            result[i]=enc2dArray;
        }
        return result;
    }
}
