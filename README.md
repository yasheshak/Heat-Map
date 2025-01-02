# **Gene Phase and Expression Heatmap Generator**

## **Overview**
This tool generates a heatmap visualizing the correlation between gene peak phases (CT values) and their normalized expression levels. Using input files containing phase data and expression data, the script produces a graphical output to explore patterns in gene expression.

---

## **Features**
1. **Data Input**:
   - Reads phase data (`.phase`) containing gene IDs and peak phases (CT values).
   - Reads expression data (`.exp`) with gene IDs and FPKM values across multiple conditions.

2. **Data Normalization**:
   - Normalizes gene expression data for comparability.
   - Maps normalized values to a visually distinct color scale.

3. **Heatmap Generation**:
   - Produces a heatmap representing gene expression levels sorted by peak phase (CT values).
   - Includes a color bar for interpreting normalized expression values.

4. **Customizable Output**:
   - Saves the heatmap as a high-resolution PNG file.

---

## **Dependencies**
- **Python 3.x**
- **Matplotlib** (`pip install matplotlib`)
- **NumPy** (`pip install numpy`)

---

## **Usage**
1. **Input Files**:
   - **Phase File (`.phase`)**: Gene IDs and peak phases (CT values).
   - **Expression File (`.exp`)**: Gene IDs and FPKM values across conditions.

2. **Execution**:
   Run the script using the command line:
   ```bash
   python heatmap_generator.py -p <phase_file> -c <expression_file> -o <output_file>
