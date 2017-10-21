/*
 * SaveImage.java
 *
 * Created on Nedela, 2008, marec 30, 3:08
 *
 * SK: trieda sluzi na ukladanie obrazkov (RenderedImage) do konkretneho suboru
 * En: class to save Images to files
 */

package InOut;

import java.io.FileWriter;
import java.io.*;
import java.awt.*;
import java.awt.image.*;
import java.lang.Object;
import java.io.File;
import javax.imageio.ImageIO;

/**
 *
 * @author Jiri Silha
 */
public class SaveImage {
    
    /** Creates a new instance of SaveImage */
    public SaveImage() {
    }
    
    /**
     * savePictureJpg()
     *
     * Method to save picture to JPG format.
     */
	
    public void savePictureJpg(String name, RenderedImage image ){
	try{
		File file=new File(name+".jpg");
		ImageIO.write(image, "jpg", file);
	}
	catch (IOException e){}
    }
    
    /**
     * savePictureBmp()
     *
     * Method to save picture to JPG format.
     */
	
    public void savePictureBmp(String name, RenderedImage image ){
	try{
		File file=new File(name+".bmp");
		ImageIO.write(image, "bmp", file);
	}
	catch (IOException e){}
    }

    /**
     * savePicturePng()
     *
     * Method to save picture to JPG format.
     */

    public void savePicturePng(String name, RenderedImage image ){
	try{
		File file=new File(name+".png");
		ImageIO.write(image, "png", file);
	}
	catch (IOException e){}
    }
}
