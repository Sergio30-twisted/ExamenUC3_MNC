/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/GUIForms/JFrame.java to edit this template
 */
package pruebaregresionmincuadgui;

import javax.swing.table.DefaultTableModel;
import java.text.DecimalFormat;
import javax.swing.JOptionPane;
import javax.swing.table.DefaultTableModel;

/**
 *
 * @author Sergio Aboytia
 */
public class MainFrame extends javax.swing.JFrame {
    
    DefaultTableModel modelo;

    private static final java.util.logging.Logger logger = java.util.logging.Logger.getLogger(MainFrame.class.getName());

    /**
     * Creates new form MainFrame
     */
    public MainFrame() {
        initComponents();
        modelo = new DefaultTableModel();
        tblDatos.setModel(modelo);
    }

    private class SolveResult {
        double[] solution;
        String steps;
        SolveResult(double[] sol, String s) { solution = sol; steps = s; }
    }
    
    private SolveResult gaussJordanEliminacion(double[][] A) {
        // A es la matriz aumentada n x (n+1)
        int n = A.length;
        StringBuilder log = new StringBuilder();
        DecimalFormat df = new DecimalFormat("#.######");
        // Copiar para no modificar original si hace falta (ya trabajamos en copia)
        double[][] a = new double[n][n+1];
        for (int i = 0; i < n; i++)
            System.arraycopy(A[i], 0, a[i], 0, n + 1);
            log.append("Matriz aumentada inicial:\n");
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n + 1; j++) {
                log.append(df.format(a[i][j])).append("\t");
            }
            log.append("\n");
        }
        log.append("\n");

        // Eliminación Gauss-Jordan con pivoteo parcial
        for (int col = 0; col < n; col++) {
            // Pivoteo parcial: buscar fila con mayor |valor| en columna col (desde col hacia abajo)
            int pivotRow = col;
            double maxVal = Math.abs(a[col][col]);
            for (int r = col + 1; r < n; r++) {
                if (Math.abs(a[r][col]) > maxVal) {
                    maxVal = Math.abs(a[r][col]);
                    pivotRow = r;
                }
            }
            if (Math.abs(a[pivotRow][col]) < 1e-12) {
                log.append("Pivote prácticamente cero en columna ").append(col).append(". El sistema puede ser singular.\n");
                return new SolveResult(null, log.toString());
            }
            // Intercambiar filas si es necesario
            if (pivotRow != col) {
                double[] tmp = a[col];
                a[col] = a[pivotRow];
                a[pivotRow] = tmp;
                log.append("Intercambio fila ").append(col).append(" con fila ").append(pivotRow).append("\n");
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n + 1; j++) {
                        log.append(df.format(a[i][j])).append("\t");
                    }
                    log.append("\n");
                }
                log.append("\n");
            }

            // Normalizar fila del pivote (hacer 1 el pivote)
            double piv = a[col][col];
            log.append("Normalizando fila ").append(col).append(" (dividir por ").append(df.format(piv)).append(")\n");
            for (int j = col; j < n + 1; j++) a[col][j] /= piv;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n + 1; j++) {
                    log.append(df.format(a[i][j])).append("\t");
                }
                log.append("\n");
            }
            log.append("\n");

            // Eliminar otras filas (hacer 0 en columna col)
            for (int row = 0; row < n; row++) {
                if (row == col) continue;
                double factor = a[row][col];
                if (Math.abs(factor) < 1e-15) continue;
                log.append("Hacer fila ").append(row).append(" = fila ").append(row).append(" - (").append(df.format(factor)).append(") * fila ").append(col).append("\n");
                for (int j = col; j < n + 1; j++) {
                    a[row][j] -= factor * a[col][j];
                }
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n + 1; j++) {
                        log.append(df.format(a[i][j])).append("\t");
                    }
                    log.append("\n");
                }
                log.append("\n");
            }
        }

        // Extraer solución
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = a[i][n]; // última columna
        }
        log.append("Solución obtenida:\n");
        for (int i = 0; i < n; i++) log.append("x").append(i).append(" = ").append(df.format(x[i])).append("\n");

        return new SolveResult(x, log.toString());
}
    
    private DefaultTableModel modelo() {
        return (DefaultTableModel) tblDatos.getModel(); 
}
    
    private void calcularRegresionPolinomial() {
    try {
        DefaultTableModel m = modelo();
        int n = m.getRowCount();
        if (n < 2) { txtResultados.setText("Se requieren al menos 2 puntos."); return; }
        int grado = Integer.parseInt(txtGrado.getText().trim());
        if (grado < 1) { txtResultados.setText("El grado debe ser >= 1."); return; }
        // Usamos X en columna 0 y Y en la última
        int colX = 0;
        int colY = m.getColumnCount() - 1;

        // Precalcular potencias de X: S[k] = Σ x^k para k=0..2grado
        double[] S = new double[2 * grado + 1];
        double[] T = new double[grado + 1]; // T[k] = Σ y * x^k
        for (int i = 0; i < n; i++) {
            double x = Double.parseDouble(String.valueOf(m.getValueAt(i, colX)));
            double y = Double.parseDouble(String.valueOf(m.getValueAt(i, colY)));
            double xp = 1.0;
            for (int k = 0; k <= 2 * grado; k++) {
                S[k] += xp;
                xp *= x;
            }
            double xp2 = 1.0;
            for (int k = 0; k <= grado; k++) {
                T[k] += y * xp2;
                xp2 *= x;
            }
        }

        // Construir matriz aumentada (grado+1)x(grado+2)
        int dim = grado + 1;
        double[][] A = new double[dim][dim + 1];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                A[i][j] = S[i + j]; // Σ x^(i+j)
            }
            A[i][dim] = T[i]; // Σ y x^i
        }

        // Registrar pasos parciales (sumas) antes de resolver
        StringBuilder pre = new StringBuilder();
        DecimalFormat df = new DecimalFormat("#.######");
        pre.append("REGRESIÓN POLINOMIAL grado ").append(grado).append("\n\n");
        pre.append("Sumas S_k = Σ x^k:\n");
        for (int k = 0; k < S.length; k++) pre.append("S[").append(k).append("] = ").append(df.format(S[k])).append("\n");
        pre.append("\nT_k = Σ y x^k:\n");
        for (int k = 0; k < T.length; k++) pre.append("T[").append(k).append("] = ").append(df.format(T[k])).append("\n\n");
        pre.append("Matriz aumentada de las ecuaciones normales:\n");
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim + 1; j++) pre.append(df.format(A[i][j])).append("\t");
            pre.append("\n");
        }
        pre.append("\n");

        SolveResult res = gaussJordanEliminacion(A);
        if (res.solution == null) {
            txtResultados.setText(pre.toString() + res.steps);
            return;
        }

        // Formatear ecuación polinomial
        StringBuilder ecu = new StringBuilder();
        DecimalFormat df2 = new DecimalFormat("#.######");
        ecu.append("Polinomio resultante:\n");
        for (int i = 0; i < res.solution.length; i++) {
            double coef = res.solution[i];
            if (i == 0) ecu.append(df2.format(coef));
            else {
                ecu.append(" + ").append(df2.format(coef)).append(" x");
                if (i > 1) ecu.append("^").append(i);
            }
        }
        // Resultado final
        txtResultados.setText(pre.toString() + res.steps + "\n" + ecu.toString());
    } catch (Exception ex) {
        txtResultados.setText("Error al calcular regresión polinomial: " + ex.getMessage());
    }
}
    
    private void calcularRegresionMultiple() {
    try {
        DefaultTableModel m = modelo();
        int n = m.getRowCount();
        int columnas = m.getColumnCount();
        if (n < 2) { txtResultados.setText("Se requieren al menos 2 puntos."); return; }
        // Última columna es Y; las anteriores (0.. columnas-2) son X1..Xm
        int mVars = columnas - 1; // número de variables independientes
        if (mVars < 1) { txtResultados.setText("No hay variables independientes."); return; }

        // Construir XᴛX (dim x dim) y Xᴛy (dim)
        // dim = mVars + 1 (intercepto)
        int dim = mVars + 1;
        double[][] normal = new double[dim][dim + 1]; // aumentada

        // Llenar normal: normal[i][j] = Σ (x_i * x_j)
        // donde x_0 = 1 (intercepto), x_k = valor de variable k (k=1..mVars)
        for (int i = 0; i < n; i++) {
            double[] rowX = new double[dim];
            rowX[0] = 1.0;
            for (int j = 1; j < dim; j++) {
                rowX[j] = Double.parseDouble(String.valueOf(m.getValueAt(i, j - 1)));
            }
            double y = Double.parseDouble(String.valueOf(m.getValueAt(i, columnas - 1)));
            for (int r = 0; r < dim; r++) {
                for (int c = 0; c < dim; c++) {
                    normal[r][c] += rowX[r] * rowX[c];
                }
                normal[r][dim] += rowX[r] * y; // columna aumentada
            }
        }

        // Mostrar la matriz normal antes de resolver
        StringBuilder pre = new StringBuilder();
        DecimalFormat df = new DecimalFormat("#.######");
        pre.append("REGRESIÓN LINEAL MÚLTIPLE\n");
        pre.append("Variables independientes (m) = ").append(mVars).append("\n\n");
        pre.append("Matriz aumentada X^T * X | X^T * y :\n");
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim + 1; j++) {
                pre.append(df.format(normal[i][j])).append("\t");
            }
            pre.append("\n");
        }
        pre.append("\n");

        SolveResult res = gaussJordanEliminacion(normal);
        if (res.solution == null) { txtResultados.setText(pre.toString() + res.steps); return; }

        // Formatear ecuación Y = b0 + b1 x1 + b2 x2 + ...
        StringBuilder ecu = new StringBuilder();
        DecimalFormat df2 = new DecimalFormat("#.######");
        ecu.append("Ecuación resultante:\nY = ").append(df2.format(res.solution[0]));
        for (int i = 1; i < res.solution.length; i++) {
            ecu.append(" + ").append(df2.format(res.solution[i])).append(" X").append(i);
        }

        txtResultados.setText(pre.toString() + res.steps + "\n" + ecu.toString());
    } catch (Exception ex) {
        txtResultados.setText("Error al calcular regresión múltiple: " + ex.getMessage());
    }
}
    
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new javax.swing.JLabel();
        cmbMetodo = new javax.swing.JComboBox<>();
        jLabel2 = new javax.swing.JLabel();
        txtPuntos = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        txtVariables = new javax.swing.JTextField();
        jScrollPane1 = new javax.swing.JScrollPane();
        tblDatos = new javax.swing.JTable();
        btnCalcular = new javax.swing.JButton();
        btnLimpiar = new javax.swing.JButton();
        btnSalir = new javax.swing.JButton();
        jLabel4 = new javax.swing.JLabel();
        jScrollPane2 = new javax.swing.JScrollPane();
        txtResultados = new javax.swing.JTextArea();
        jLabel5 = new javax.swing.JLabel();
        txtGrado = new javax.swing.JTextField();
        btnGenerar = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setResizable(false);

        jLabel1.setText("Seleccione el método a emplear:");

        cmbMetodo.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Regresión Lineal Simple", "Regresión Polinomial", "Regresión Lineal Múltiple", " " }));

        jLabel2.setText("Cantidad de puntos:");

        txtPuntos.setText("00000000");

        jLabel3.setText("Cantidad de variables independientes:");

        txtVariables.setText("00000000");
        txtVariables.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtVariablesActionPerformed(evt);
            }
        });

        tblDatos.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                {null, null, null, null},
                {null, null, null, null},
                {null, null, null, null},
                {null, null, null, null}
            },
            new String [] {
                "Title 1", "Title 2", "Title 3", "Title 4"
            }
        ) {
            boolean[] canEdit = new boolean [] {
                false, false, false, false
            };

            public boolean isCellEditable(int rowIndex, int columnIndex) {
                return canEdit [columnIndex];
            }
        });
        jScrollPane1.setViewportView(tblDatos);

        btnCalcular.setIcon(new javax.swing.ImageIcon(getClass().getResource("/resources/calculatepin.png"))); // NOI18N
        btnCalcular.setText("Calcular");
        btnCalcular.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcularActionPerformed(evt);
            }
        });

        btnLimpiar.setIcon(new javax.swing.ImageIcon(getClass().getResource("/resources/clearpin.png"))); // NOI18N
        btnLimpiar.setText("Limpiar");
        btnLimpiar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnLimpiarActionPerformed(evt);
            }
        });

        btnSalir.setIcon(new javax.swing.ImageIcon(getClass().getResource("/resources/exitpin.png"))); // NOI18N
        btnSalir.setText("Salir");
        btnSalir.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSalirActionPerformed(evt);
            }
        });

        jLabel4.setText("Resultado");

        txtResultados.setColumns(20);
        txtResultados.setRows(5);
        jScrollPane2.setViewportView(txtResultados);

        jLabel5.setText("Grado del polinomio:");

        txtGrado.setText("00000000");

        btnGenerar.setText("Generar Tabla");
        btnGenerar.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnGenerarActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addGap(18, 18, 18)
                        .addComponent(cmbMetodo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(btnGenerar)
                        .addGap(88, 88, 88))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(jLabel2)
                        .addGap(18, 18, 18)
                        .addComponent(txtPuntos, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(18, 18, 18)
                        .addComponent(jLabel3)
                        .addGap(18, 18, 18)
                        .addComponent(txtVariables, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(0, 34, Short.MAX_VALUE))))
            .addGroup(layout.createSequentialGroup()
                .addGap(35, 35, 35)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel5)
                        .addGap(18, 18, 18)
                        .addComponent(txtGrado, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(jScrollPane1)
                        .addComponent(jScrollPane2)
                        .addGroup(layout.createSequentialGroup()
                            .addComponent(btnCalcular)
                            .addGap(87, 87, 87)
                            .addComponent(btnLimpiar)
                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(btnSalir))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jLabel4)
                .addGap(232, 232, 232))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(cmbMetodo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(txtPuntos, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3)
                    .addComponent(txtVariables, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(26, 26, 26)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel5)
                    .addComponent(txtGrado, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnGenerar))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 35, Short.MAX_VALUE)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 150, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnCalcular)
                    .addComponent(btnLimpiar)
                    .addComponent(btnSalir))
                .addGap(18, 18, 18)
                .addComponent(jLabel4)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 129, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(60, 60, 60))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnLimpiarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnLimpiarActionPerformed
        modelo.setRowCount(0);
        txtPuntos.setText("");
        txtVariables.setText("");
        txtGrado.setText("");
        txtResultados.setText("");
    }//GEN-LAST:event_btnLimpiarActionPerformed

    private void btnSalirActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSalirActionPerformed
        System.exit(0);
    }//GEN-LAST:event_btnSalirActionPerformed

    private void btnCalcularActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCalcularActionPerformed
        if (tblDatos.isEditing()) {
    tblDatos.getCellEditor().stopCellEditing();
}
        String metodo = cmbMetodo.getSelectedItem().toString();

        DefaultTableModel model = (DefaultTableModel) tblDatos.getModel();

if (model.getRowCount() == 0) {
    JOptionPane.showMessageDialog(this,
            "Primero debes presionar *Generar* para crear la tabla con los datos.",
            "Error",
            JOptionPane.ERROR_MESSAGE
    );
    return;
}
        switch (metodo) {
            case "Regresión Lineal Simple":
                calcularRegresionSimple();
                break;

            case "Regresión Polinomial":
                calcularRegresionPolinomial();
                break;

            case "Regresión Lineal Múltiple":
                calcularRegresionMultiple();
                break;
        }
    }//GEN-LAST:event_btnCalcularActionPerformed

    private void btnGenerarActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnGenerarActionPerformed
        
    try {
        int puntos = Integer.parseInt(txtPuntos.getText());
        int variables = Integer.parseInt(txtVariables.getText());

        if (puntos <= 0 || variables <= 0) {
            JOptionPane.showMessageDialog(this, "Debe ingresar valores mayores a 0.");
            return;
        }

        // Crear modelo
        DefaultTableModel modelo = new DefaultTableModel();
        tblDatos.setModel(modelo);

        // Agregar columnas X1...Xn
        for (int i = 1; i <= variables; i++) {
            modelo.addColumn("X" + i);
        }

        // Agregar columna Y (dependiente)
        modelo.addColumn("Y");

        // Agregar filas vacías según cantidad de puntos
        for (int i = 0; i < puntos; i++) {
            modelo.addRow(new Object[variables + 1]);
        }

        txtResultados.setText("✅ Tabla generada correctamente.\nIngrese los valores y luego presione CALCULAR.");
    
    } catch (NumberFormatException ex) {
        JOptionPane.showMessageDialog(this, "Debe ingresar números válidos.");
    }
    }//GEN-LAST:event_btnGenerarActionPerformed

    private void txtVariablesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtVariablesActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtVariablesActionPerformed

    private void calcularRegresionSimple() {
        
    DefaultTableModel model = (DefaultTableModel) tblDatos.getModel();

    // Forzar guardar la edición de la celda activa
    if (tblDatos.isEditing()) {
        tblDatos.getCellEditor().stopCellEditing();
    }

    // Validar que existan filas
    if (model.getRowCount() == 0) {
        JOptionPane.showMessageDialog(this,
            "No hay datos en la tabla.\nPrimero presiona GENERAR y luego llena la tabla.",
            "Error",
            JOptionPane.ERROR_MESSAGE);
        return;
    }

    int n = model.getRowCount();
    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;

    for (int i = 0; i < n; i++) {

        // Validar celdas vacías antes de leer
        if (model.getValueAt(i, 0) == null || model.getValueAt(i, 1) == null) {
            JOptionPane.showMessageDialog(this,
                "Hay celdas vacías en la fila " + (i + 1),
                "Datos incompletos",
                JOptionPane.WARNING_MESSAGE);
            return;
        }

        double x = Double.parseDouble(model.getValueAt(i, 0).toString());
        double y = Double.parseDouble(model.getValueAt(i, 1).toString());

        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += x * x;
    }

    double b = (n * sumXY - sumX * sumY) / (n * sumX2 - Math.pow(sumX, 2));
    double a = (sumY - b * sumX) / n;

    txtResultados.setText(
        "REGRESIÓN LINEAL SIMPLE\n\n" +
        "Ecuación resultante:\n" +
        "Y = " + a + " + " + b + "X\n\n"
    );
}

    
    private void txtVariablesFocusLost(java.awt.event.FocusEvent evt) {                                     
    int variables = Integer.parseInt(txtVariables.getText());
    modelo.setColumnCount(0);  // Limpiar columnas previas

    // Generar columnas dinámicamente
    for (int i = 1; i <= variables; i++) {
        modelo.addColumn("X" + i);
    }

    modelo.addColumn("Y"); // Última columna siempre Y (dependiente)
}
    
    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ReflectiveOperationException | javax.swing.UnsupportedLookAndFeelException ex) {
            logger.log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(() -> new MainFrame().setVisible(true));
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnCalcular;
    private javax.swing.JButton btnGenerar;
    private javax.swing.JButton btnLimpiar;
    private javax.swing.JButton btnSalir;
    private javax.swing.JComboBox<String> cmbMetodo;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JTable tblDatos;
    private javax.swing.JTextField txtGrado;
    private javax.swing.JTextField txtPuntos;
    private javax.swing.JTextArea txtResultados;
    private javax.swing.JTextField txtVariables;
    // End of variables declaration//GEN-END:variables
}
