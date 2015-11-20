/**

Title:       PlotReader.java
Author:      Nicholas Geary (adapted from MainClass.java by Brian Coopersmith)
Description: Reads data points from an image file and saves them to a text file.
To use:      1) compile with "javac PlotReader.java"
             2) run with "java PlotReader [-option] arg [-option] arg"

**Option**   **Description**                         **Input type**     **Default**
 -i           name of image file                      String             (no default)
 -o           name of save file                       String             PlotReaderOutput.txt
 -s           y-axis scale (linear or log10)          String             linear
 -xs          x-axis scale (linear or log10)          String             linear
 -nf          # format (# of decimal places in data)  int                5

**/

import java.awt.*;
import java.awt.event.*;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.swing.*;
import javax.swing.text.*;
import java.io.*;
import javax.imageio.*;
import java.text.*;

public class PlotReader extends JFrame implements MouseListener, KeyListener, ActionListener {
    
    private static File imageFile;
    private static String imageFileName;
    private static boolean isLinear = true;
    private static boolean XisLinear = true;
    private static String saveFileName;
    private static File saveFile;
    private static int numPlaces = 5;
    private String nextInstruction;
    private StringBuilder sb = new StringBuilder("Click on the low X point of the graph.");
    private int curStep = 0;
    private double lowX_pix, highX_pix, lowY_pix, highY_pix;
    private double lowX_val, highX_val, lowY_val, highY_val;
    private double[] point = new double[3];
    private ArrayList<double[]> data_points=new ArrayList<double[]>();
    private JPanel picturePanel;
    private JPanel textPanel;
    private JLabel thePicture;
    private JTextField textBox;
    private JLabel instructionsText;
    private JPanel innerPanel;
    private JLabel modeLabel;
    private JButton saveButton;
    private JPanel savePanel;
    private JLabel saveLabel;
    private static DecimalFormat df;
    private static PrintWriter pw;
    
    public static void main (String [] args) {
	boolean imageFileDefined = false;
	boolean saveFileDefined = false;
	for (int i=0; i<args.length; i++) {
	    if (args[i].toLowerCase().equals("-i")) {
		imageFileName = args[i+1];
		imageFileDefined = true;
	    }
	    if (args[i].toLowerCase().equals("-s")) {
		isLinear = !(args[i+1].toLowerCase().contains("log"));
	    }
	    if (args[i].toLowerCase().equals("-xs")) {
		XisLinear = !(args[i+1].toLowerCase().contains("log"));
	    }
	    if (args[i].toLowerCase().equals("-o")) {
		saveFileName = args[i+1];
		saveFileDefined = true;
	    }
	    if (args[i].toLowerCase().equals("-nf")) {
		numPlaces = Integer.parseInt(args[i+1]);
	    }
	}
	if (!imageFileDefined) {
	    System.out.println("Please enter an image file name. See source code of PlotReader.java for help.");
	    System.exit(1);
	}
	if (!saveFileDefined) {
	    saveFileName = new String("PlotReaderOutput.txt");
	}
	StringBuilder dfsb = new StringBuilder("#.");
	for (int i=0; i<numPlaces; i++) {
	    dfsb.append("#");
	}
	df = new DecimalFormat(dfsb.toString());
	PlotReader PR = new PlotReader();
	
    }
    
    
    public PlotReader() {
	
	imageFile = new File(imageFileName);
	
	textBox = new JTextField("");
	instructionsText = new JLabel(sb.toString());
	textPanel = new JPanel(new GridLayout(2, 1));
	instructionsText.setFont(new Font(instructionsText.getFont().getName(), Font.BOLD, 14));
	instructionsText.setForeground(Color.RED);
	textPanel.add(instructionsText);
	innerPanel = new JPanel(new GridLayout(1, 2));
	modeLabel = new JLabel();
	if (isLinear) {
	    if (XisLinear) {
		modeLabel.setText("Y-axis: linear  X-axis: linear");
	    } else {
		modeLabel.setText("Y-axis: linear  X-axis: logarithmic");
	    }
	} else {
	    if (XisLinear) {
		modeLabel.setText("Y-axis: logarithmic  X-axis: linear");
	    } else {
		modeLabel.setText("Y-axis: logarithmic  X-axis: logarithmic");
	    } 
	}
	innerPanel.add(textBox);
	innerPanel.add(modeLabel);
	textPanel.add(innerPanel);
	
	picturePanel = new JPanel (new GridLayout(1, 1));
	BufferedImage bi = null;
	try {
	    bi = ImageIO.read(imageFile);
	} catch (IOException e) {
	    e.printStackTrace();
	}
	thePicture = new JLabel(new ImageIcon(bi));
	picturePanel.add(thePicture);
	
	saveButton = new JButton("Save to file");
	saveLabel = new JLabel("");
	savePanel = new JPanel(new GridLayout(1,2));
	savePanel.add(saveButton);
	savePanel.add(saveLabel);
	
	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	setLayout(new BorderLayout());
	add(textPanel, BorderLayout.NORTH);
	add(picturePanel, BorderLayout.CENTER);
	add(savePanel, BorderLayout.SOUTH);
	pack();
	setVisible(true);
	
	this.addMouseListener(this);
	textBox.addKeyListener(this);
	saveButton.addActionListener(this);
	
    }
    
    
    public Image getImage(File imageFile) {
	BufferedImage img = null;
	try {
	    img = ImageIO.read(imageFile);
	} catch (IOException e) {
	    e.printStackTrace();
	}
	return img;
    }
    
    private void refreshText() {
	instructionsText.setText(sb.toString());
    }
    
    
    
    public void mouseClicked(MouseEvent me) {
	
	switch(curStep) {
	case 0:
	    lowX_pix = me.getX();
	    sb.setLength(0);
	    sb.append("Enter the low X in the text field below and press enter.");
	    curStep++;
	    break;
	case 2:
	    highX_pix = me.getX();
	    sb.setLength(0);
	    sb.append("Enter the high X in the text field below and press enter.");
	    curStep++;
	    break;
	case 4:
	    lowY_pix = me.getY();
	    sb.setLength(0);
	    sb.append("Enter the low Y in the text field below and press enter.");
	    curStep++;
	    break;
	case 6:
	    highY_pix = me.getY();
	    sb.setLength(0);
	    sb.append("Enter the high Y in the text field below and press enter.");
	    curStep++;
	    break;
	case 8:
	    if (XisLinear) {
		point[0] = (me.getX()-lowX_pix)*((highX_val-lowX_val)/(highX_pix-lowX_pix))+lowX_val;
	    } else {
		point[0] = Math.pow(10,(me.getX()-lowX_pix)*(Math.log10(highX_val/lowX_val)/(highX_pix-lowX_pix))+Math.log10(lowX_val));
	    }
	    point[1] = getYVal(me.getY());
	    sb.setLength(0);
	    sb.append("Click on the point of the upper error bar for the same data point.");
	    curStep++;
	    break;
	case 9:
	    point[2] = getYVal(me.getY());
	    sb.setLength(0);
	    sb.append("Click on the point of the lower error bar for a new data point.");
	    data_points.add(point);
	    point = new double[3];
	    curStep--;
	    for(int i=0;i<data_points.size();i++) {
		System.out.println(df.format(data_points.get(i)[0])+" "+df.format(data_points.get(i)[1])+" "+df.format(data_points.get(i)[2]));
	    }
	    System.out.println("");
	    break;
	}
	
	repaint();
	refreshText();
	
    }
    
    
    public void mousePressed(MouseEvent me) {}
    public void mouseReleased(MouseEvent me) {}
    public void mouseEntered(MouseEvent me) {}
    public void mouseExited(MouseEvent me) {}
    
    
    public void keyPressed(KeyEvent ke) {
	if(ke.getKeyChar()=='\n') {
	    switch(curStep)
		{
		case 1:
		    lowX_val = Double.parseDouble(textBox.getText());
		    textBox.setText("");
		    sb.setLength(0);
		    sb.append("Click on the high X point of the graph.");
		    curStep++;
		    break;
		case 3:
		    highX_val = Double.parseDouble(textBox.getText());
		    textBox.setText("");
		    sb.setLength(0);
		    sb.append("Click on the low Y point of the graph.");
		    curStep++;
		    break;
		case 5:
		    lowY_val = Double.parseDouble(textBox.getText());
		    textBox.setText("");
		    sb.setLength(0);
		    sb.append("Click on the high Y point of the graph.");
		    curStep++;
		    break;
		case 7:
		    highY_val = Double.parseDouble(textBox.getText());
		    textBox.setText("");
		    sb.setLength(0);
		    sb.append("Click on the point of the lower error bar for a data point");
		    curStep++;
		    break;
		}
	    repaint();
	    refreshText();
	}	
    }
    
    public void keyReleased(KeyEvent ke){}
    public void keyTyped(KeyEvent ke){}
    
    public void actionPerformed(ActionEvent e) {
	try {
	    saveFile = new File(saveFileName);
	    pw = new PrintWriter(saveFile);
	    pw.println("     X     |    LoY    |    HiY    ");
	    for(int i=0;i<data_points.size();i++) {
		pw.println(df.format(data_points.get(i)[0])+" "+df.format(data_points.get(i)[1])+" "+df.format(data_points.get(i)[2]));
	    }
	    pw.close();
	    saveLabel.setText("Saved data to " + saveFileName); 
	} catch (IOException ioe) {
	    System.out.println("Could not save to file.");
	    saveLabel.setText("Could not save to file.");
	}
    }
    
    private double getYVal(double pixel) {
	if (isLinear) {
	    return ((pixel-lowY_pix)*((highY_val-lowY_val)/(highY_pix-lowY_pix))+lowY_val);
	} else {
	    double temp = (pixel-lowY_pix)*(Math.log10(highY_val/lowY_val)/(highY_pix-lowY_pix))+Math.log10(lowY_val);
	    return Math.pow(10, temp);
	}
    }
    
    
}