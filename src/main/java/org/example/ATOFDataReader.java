package org.example;

import org.jlab.jnp.hipo4.data.*;
import org.jlab.jnp.hipo4.io.HipoReader;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.util.*;

public class ATOFDataReader {

    private static final float VEFF = 20.0f; // Speed of light in cm/ns
    private static final float WEDGE_SPACING = 3.0f; // Wedge spacing in mm
    private static final float Z_TOLERANCE = 5.0f; // Tolerance for Z matching
    private static final int WEDGES_PER_BAR = 10;
    private static final int NUM_BARS = 60;

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println("Usage: java ATOFDataReader <input.hipo>");
            System.exit(1);
        }

        String inputHipoFile = args[0];
        HipoReader reader = new HipoReader();
        reader.open(inputHipoFile);

        SchemaFactory schemaFactory = reader.getSchemaFactory();
        Bank adcBank = new Bank(schemaFactory.getSchema("ATOF::adc"));
        Bank mcTrueBank = new Bank(schemaFactory.getSchema("MC::True"));
        Event event = new Event();


        XYSeries zTruthVsZBar = new XYSeries("ZTruth vs ZBar");
        XYSeries zTruthVsZWedge = new XYSeries("ZTruth vs ZWedge");
        XYSeries zBarVsZWedge = new XYSeries("ZBar vs ZWedge");
        XYSeries deltaZBar = new XYSeries("Delta Z (Truth - ZBar)");
        XYSeries deltaZWedge = new XYSeries("Delta Z (Truth - ZWedge)");

        int eventId = 0;

        while (reader.hasNext()) {
            reader.nextEvent(event);
            event.read(adcBank);
            event.read(mcTrueBank);
            eventId++;

            List<Hit> truthBarHits = new ArrayList<>();
            List<Hit> truthWedgeHits = new ArrayList<>();
            extractTruthHits(mcTrueBank, truthBarHits, truthWedgeHits);

            List<Hit> barHits = new ArrayList<>();
            List<Hit> wedgeHits = new ArrayList<>();
            extractADCData(adcBank, barHits, wedgeHits);


            List<Hit> zBarHits = calculateZBar(barHits);

	    for (Hit truthBar : truthBarHits) {
                for (Hit bar : zBarHits) {
                    zTruthVsZBar.add(truthBar.z, bar.z);
                    deltaZBar.add(eventId, truthBar.z - bar.z);
                }
            }

            for (Hit truthWedge : truthWedgeHits) {
                for (Hit wedge : wedgeHits) {
                    if (Math.abs(truthWedge.z - wedge.z) < Z_TOLERANCE) {
                        zTruthVsZWedge.add(truthWedge.z, wedge.z);
                        deltaZWedge.add(eventId, truthWedge.z - wedge.z);
                    }
                }
            }
            for (Hit bar : zBarHits) {
                for (Hit wedge : wedgeHits) {
                    if (Math.abs(bar.z - wedge.z) < Z_TOLERANCE) {
                        zBarVsZWedge.add(bar.z, wedge.z);
                    }
                }
            }
            printClusters(eventId, zBarHits, wedgeHits);
        }

        reader.close();

        SwingUtilities.invokeLater(() -> {
            showScatterPlot(zTruthVsZBar, "ZTruth vs ZBar", "ZTruth", "ZBar");
            showScatterPlot(zTruthVsZWedge, "ZTruth vs ZWedge", "ZTruth", "ZWedge");
            showScatterPlot(zBarVsZWedge, "ZBar vs ZWedge", "ZBar", "ZWedge");
            showScatterPlot(deltaZBar, "Delta Z (Truth - ZBar)", "Event ID", "Delta Z");
            showScatterPlot(deltaZWedge, "Delta Z (Truth - ZWedge)", "Event ID", "Delta Z");
        });

        System.out.println("Processing Complete.");
    }

    private static void extractTruthHits(Bank mcTrueBank, List<Hit> barTruthHits, List<Hit> wedgeTruthHits) {
        for (int i = 0; i < mcTrueBank.getRows(); i++) {
            int detector = mcTrueBank.getInt("detector", i);
            int hitn = mcTrueBank.getInt("hitn", i);
            float avgZ = mcTrueBank.getFloat("avgZ", i);
            float avgT = mcTrueBank.getFloat("avgT", i);
            float px = mcTrueBank.getFloat("px", i);
            float py = mcTrueBank.getFloat("py", i);
            float phi = (float) Math.atan2(py, px);

            if (detector == 24) {
                if (hitn == 1) barTruthHits.add(new Hit(avgZ, phi, avgT));
                if (hitn == 2) wedgeTruthHits.add(new Hit(avgZ, phi, avgT));
            }
        }
    }

    private static void extractADCData(Bank adcBank, List<Hit> barHits, List<Hit> wedgeHits) {
        for (int i = 0; i < adcBank.getRows(); i++) {
            int layer = adcBank.getByte("layer", i);
            int component = adcBank.getShort("component", i);
            int order = adcBank.getByte("order", i);
            float time = adcBank.getFloat("time", i);

            if (layer == 0) barHits.add(new Hit(component, order, time));
            if (layer >= 10 && layer <= 19) {
                float zWedge = calculateZWedge(component);
                float phiWedge = calculatePhi(component);
                wedgeHits.add(new Hit(zWedge, phiWedge, time));
            }
        }
    }

    private static List<Hit> calculateZBar(List<Hit> barHits) {
        Map<Integer, Float> tLeft = new HashMap<>();
        Map<Integer, Float> tRight = new HashMap<>();
        List<Hit> zBarHits = new ArrayList<>();

        for (Hit hit : barHits) {
            if (hit.order == 0) tLeft.put(hit.component, hit.time);
            if (hit.order == 1) tRight.put(hit.component, hit.time);
        }

        for (int component : tLeft.keySet()) {
            if (tRight.containsKey(component)) {
                float zBar = VEFF * (tRight.get(component) - tLeft.get(component)) / 2;
                float phiBar = calculatePhi(component);
                float avgTime = (tLeft.get(component) + tRight.get(component)) / 2;
                zBarHits.add(new Hit(zBar, phiBar, avgTime));
            }
        }
        return zBarHits;
    }

    private static float calculateZWedge(int component) {
        return (component % WEDGES_PER_BAR - (WEDGES_PER_BAR - 1) / 2.0f) * WEDGE_SPACING;
    }

    private static float calculatePhi(int component) {
        return (float) (2.0 * Math.PI * component / NUM_BARS);
    }

    private static void printClusters(int eventId, List<Hit> bars, List<Hit> wedges) {
        System.out.printf("Event %d:\n", eventId);
        for (Hit bar : bars) System.out.printf("Bar: Z=%.2f, Phi=%.4f, Time=%.2f\n", bar.z, bar.phi, bar.time);
        for (Hit wedge : wedges) System.out.printf("Wedge: Z=%.2f, Phi=%.4f, Time=%.2f\n", wedge.z, wedge.phi, wedge.time);
    }

    private static void showScatterPlot(XYSeries series, String title, String xLabel, String yLabel) {
        XYSeriesCollection dataset = new XYSeriesCollection(series);
        JFreeChart chart = ChartFactory.createScatterPlot(title, xLabel, yLabel, dataset);
        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(new ChartPanel(chart));
        frame.pack();
        frame.setVisible(true);
    }

    static class Hit {
        float z, phi, time;
        int component, order;

        Hit(float z, float phi, float time) {
            this.z = z;
            this.phi = phi;
            this.time = time;
        }

        Hit(int component, int order, float time) {
            this.component = component;
            this.order = order;
            this.time = time;
        }
    }
}

