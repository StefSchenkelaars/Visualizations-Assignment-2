/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.List;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.Volume;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;

    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        volume = vol;

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        panel.setTransferFunctionEditor(tfEditor);


    }

    @Override
    public void changed() {
        for (TFChangeListener listener : listeners) {
            listener.changed();
        }
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    // get a voxel from the volume data by nearest neighbor interpolation
    short getVoxel(double[] coord) {

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);

        if ((x >= 0) && (x < volume.getDimX()) && (y >= 0) && (y < volume.getDimY())
                && (z >= 0) && (z < volume.getDimZ())) {
            return volume.getVoxel(x, y, z);
        } else {
            return 0;
        }
    }
    
    short getTriLinear(double[] coord) {
        int dim = 1;
        
        // Check if entire volume (dim * dim) lies inside volume
        if (( coord[0] >= dim) && (coord[0] + dim < volume.getDimX()) &&
            ( coord[1] >= dim) && (coord[1] + dim < volume.getDimY()) &&
            ( coord[2] >= dim) && (coord[2] + dim < volume.getDimZ())){
             
            // This method is from wikipedia Trilinear_interpolation
            // Calculate lattice points
            int x0 = (int) Math.floor(coord[0]);
            int x1 = x0 + 1;
            int y0 = (int) Math.floor(coord[1]);
            int y1 = y0 + 1;
            int z0 = (int) Math.floor(coord[2]);
            int z1 = z0 + 1;
            
            // Calculate distance to lattice points
            double xd = (coord[0] - x0)/(x1 - x0);
            double yd = (coord[1] - y0)/(y1 - y0);
            double zd = (coord[2] - z0)/(z1 - z0);
            
            // Interpolate along x
            double c00 = volume.getVoxel(x0, y0, z0) * (1 - xd) + volume.getVoxel(x1, y0, z0) * xd;
            double c10 = volume.getVoxel(x0, y1, z0) * (1 - xd) + volume.getVoxel(x1, y1, z0) * xd;
            double c01 = volume.getVoxel(x0, y0, z1) * (1 - xd) + volume.getVoxel(x1, y0, z1) * xd;
            double c11 = volume.getVoxel(x0, y1, z1) * (1 - xd) + volume.getVoxel(x1, y1, z1) * xd;
            
            // Interpolate along y
            double c0 = c00 * (1 - yd) + c10 * yd;
            double c1 = c01 * (1 - yd) + c11 * yd;
            
            // Interpolate along z and return
            return (short) (c0 * (1 - zd) + c1 * zd);
            
        } else {
            return 0;
        }
            
        
    }
    
    private void clearImage() {
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
    }
    
    private void setPixel(int x, int y, int pixelColor, int resolution){
        for (int i = x; i < x + resolution; i++){
            for (int j = y; j < y + resolution; j++){
                // Check if pixel still inside image, prevent out of bounds
                if (i < image.getWidth() && j < image.getHeight() ){
                    image.setRGB(i, j, pixelColor);
                }
            }
        }
    }
    
    private int toARGB(int alpha, int red, int green, int blue){
        return (alpha << 24) | (red << 16) | (green << 8) | blue;
    }
    
    void mip(double[] viewMatrix, int resolution) {

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        for (int j = 0; j < image.getHeight(); j += resolution) {
            for (int i = 0; i < image.getWidth(); i += resolution) {
                
                // Default value is 0
                int val = 0;

                // Loop through the viewVector
                // First calculate where the volume stops
                int maxZIndex = Integer.MAX_VALUE;
                for(int xyz = 0; xyz <= 2; xyz++){
                    if (viewVec[xyz] != 0.0){
                        int temp = Math.abs((int) Math.floor((uVec[xyz] * (imageCenter - i) + vVec[xyz] * (imageCenter - j) - volumeCenter[2]) / viewVec[xyz]));
                        if (temp < maxZIndex){
                            maxZIndex = temp;
                        }
                    }
                }
                
                for (int k = - maxZIndex; k < maxZIndex; k++) {
                    for(int xyz = 0; xyz <= 2; xyz++){
                        pixelCoord[xyz] = uVec[xyz] * (i - imageCenter) + vVec[xyz] * (j - imageCenter) + viewVec[xyz] * k + volumeCenter[xyz];     
                    }
                    
                    int tempVal = getTriLinear(pixelCoord);
                    if (tempVal > val){
                        val = tempVal;
                    }
                }
                
                // Normalize to max
                val = (int) Math.floor(val/max * 255);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int pixelColor = toARGB(255, val, val, val);
                setPixel(i,j,pixelColor,resolution);
            }
        }


    }
    
    void compositing(double[] viewMatrix, int resolution) {
        
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        for (int j = 0; j < image.getHeight(); j += resolution) {
            for (int i = 0; i < image.getWidth(); i += resolution) {
                
                // Calculate where the volume stops
                int maxZIndex = Integer.MAX_VALUE;
                for(int xyz = 0; xyz <= 2; xyz++) {
                    if (viewVec[xyz] != 0.0) {
                        int temp = Math.abs((int) Math.floor((uVec[xyz] * (imageCenter - i) + vVec[xyz] * (imageCenter - j) - volumeCenter[2]) / viewVec[xyz]));
                        if (temp < maxZIndex) {
                            maxZIndex = temp;
                        }
                    }
                }
                
                // Initialize color array
                ArrayList<TFColor> colors = new ArrayList<TFColor>();
                
                // Loop through the viewVector
                for (int k = - maxZIndex; k < maxZIndex; k += 10) {
                    for(int xyz = 0; xyz <= 2; xyz++) {
                        pixelCoord[xyz] = uVec[xyz] * (i - imageCenter) + vVec[xyz] * (j - imageCenter) + viewVec[xyz] * k + volumeCenter[xyz];     
                    }
                    
                    int tempVal = getVoxel(pixelCoord);
                    TFColor tempColor = tFunc.getColor(tempVal);
                    colors.add(tempColor);
                }
                
                // Calculate mean of colors array
                double d_alpha = 0.0, d_red = 0.0, d_green = 0.0, d_blue = 0.0;
                if (colors.size() > 0) {
                    ArrayList<Double> alphas = new ArrayList<Double>();
                    alphas.add(1.0 - colors.get(0).a);
                    for (int l = 1; l < colors.size(); l++) {
                        alphas.add((1.0 - colors.get(l).a) * alphas.get(l - 1));
                    }
                    d_alpha = 1.0 - alphas.get(alphas.size() - 1);

                    d_red = colors.get(0).r;
                    d_green = colors.get(0).g;
                    d_blue = colors.get(0).b;
                    for (int l = 1; l < colors.size(); l++) {
                        d_red += colors.get(l).r * alphas.get(l - 1);
                        d_green += colors.get(l).g * alphas.get(l - 1);
                        d_blue += colors.get(l).b * alphas.get(l - 1);
                    }
                }
                
                // Convert double value to color int
                int c_alpha = d_alpha <= 1.0 ? (int) Math.floor(d_alpha * 255) : 255;
                int c_red = d_red <= 1.0 ? (int) Math.floor(d_red * 255) : 255;
                int c_green = d_green <= 1.0 ? (int) Math.floor(d_green * 255) : 255;
                int c_blue = d_blue <= 1.0 ? (int) Math.floor(d_blue * 255) : 255;
                int pixelColor = toARGB(c_alpha, c_red, c_green, c_blue);
                setPixel(i,j,pixelColor,resolution);
            }
        }
        
    }
    
    void slicer(double[] viewMatrix) {

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];
                
                int val = getVoxel(pixelCoord);
                // Apply the transfer function to obtain a color
                TFColor voxelColor = tFunc.getColor(val);
                
//                System.out.println("Pixel cord: " + pixelCoord[0] + ", " + pixelCoord[1] + ", " + pixelCoord[2] + " Val: " + val);
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = toARGB(c_alpha, c_red, c_green, c_blue);
                image.setRGB(i, j, pixelColor);
            }
        }


    }

    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL2.GL_LINE_SMOOTH);
        gl.glHint(GL2.GL_LINE_SMOOTH_HINT, GL2.GL_NICEST);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL2.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL2.GL_LINE_SMOOTH);
        gl.glDisable(GL2.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {

        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        
        clearImage();
        // Resolution is the number of pixels that have the same value
        // If resolution is 3, a square of (3x3) pixels have the same value
        int resolution = 1;
        String castType = panel.getCastType();
        if (castType == "Slicer"){
            slicer(viewMatrix);
        } else if (castType == "MIP"){
            // Set resolution dependent on if currently moving
            resolution = this.interactiveMode ? 3 : 1;
            mip(viewMatrix, resolution);
        } else if (castType == "Compositing"){
            // Set resolution dependent on if currently moving
            resolution = this.interactiveMode ? 2 : 1;
            compositing(viewMatrix, resolution);
        }
        
  
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        
        // Set info on panel
        panel.setSpeedLabel(Double.toString(runningTime));
        int xResolution = image.getWidth() / resolution;
        int yResolution = image.getHeight() / resolution;
        panel.setResolutionLabel(Integer.toString(xResolution) + " x " + Integer.toString(yResolution));
        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL2.GL_BLEND);
        gl.glBlendFunc(GL2.GL_SRC_ALPHA, GL2.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
}
