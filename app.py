import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from enzyme_kinetics_fitter import fit_michaelis_menten, michaelis_menten

st.set_page_config(page_title="Catalytix: Enzyme Kinetics", page_icon="🧪")

st.title("🧪 Catalytix")
st.markdown("""
This tool fits **Michaelis–Menten** enzyme kinetics data.
Upload a CSV file with substrate concentration $[S]$ and initial velocity $v$.
""")

# Sidebar for inputs
st.sidebar.header("Upload Data")
uploaded_file = st.sidebar.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        st.write("### Data Preview")
        st.dataframe(df.head())

        # Column selection
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        
        if len(numeric_cols) < 2:
            st.error("CSV must contain at least two numeric columns.")
        else:
            st.sidebar.header("Select Columns")
            s_col = st.sidebar.selectbox("Substrate Concentration ([S])", numeric_cols, index=0)
            v_col = st.sidebar.selectbox("Velocity (v)", numeric_cols, index=1)

            if st.sidebar.button("Fit Model"):
                S = df[s_col].to_numpy(dtype=float)
                v = df[v_col].to_numpy(dtype=float)

                try:
                    fit = fit_michaelis_menten(S, v)

                    # Display Results
                    st.success("Fitting successful!")
                    
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Vmax", f"{fit.Vmax:.4g}", delta_color="off")
                    col2.metric("Km", f"{fit.Km:.4g}", delta_color="off")
                    if fit.r_squared is not None:
                        col3.metric("R²", f"{fit.r_squared:.4f}", delta_color="off")

                    # Detailed stats
                    with st.expander("See detailed statistics"):
                        st.write(f"**Vmax Std Dev:** {fit.Vmax_std:.4g}" if fit.Vmax_std else "Vmax Std Dev: N/A")
                        st.write(f"**Km Std Dev:** {fit.Km_std:.4g}" if fit.Km_std else "Km Std Dev: N/A")

                    # Plotting
                    st.write("### Fitted Curve")
                    
                    # Generate plot data
                    S_range = np.linspace(0, S.max() * 1.1, 200)
                    v_fit_curve = michaelis_menten(S_range, fit.Vmax, fit.Km)
                    
                    fig, ax = plt.subplots(figsize=(8, 5))
                    ax.scatter(S, v, label="Experimental Data", color="black", zorder=5)
                    ax.plot(S_range, v_fit_curve, label=f"Fit (Vmax={fit.Vmax:.2f}, Km={fit.Km:.2f})", color="red", linewidth=2)
                    
                    ax.set_xlabel(s_col)
                    ax.set_ylabel(v_col)
                    ax.set_title("Michaelis–Menten Fit")
                    ax.legend()
                    ax.grid(True, linestyle="--", alpha=0.5)
                    
                    st.pyplot(fig)

                    # Residuals
                    st.write("### Residuals")
                    v_pred = michaelis_menten(S, fit.Vmax, fit.Km)
                    residuals = v - v_pred
                    
                    fig_res, ax_res = plt.subplots(figsize=(8, 3))
                    ax_res.axhline(0, linestyle="--", color="gray")
                    ax_res.scatter(S, residuals, color="blue")
                    ax_res.set_xlabel(s_col)
                    ax_res.set_ylabel("Residuals")
                    ax_res.grid(True, linestyle="--", alpha=0.5)
                    
                    st.pyplot(fig_res)

                except Exception as e:
                    st.error(f"Fitting failed: {e}")

    except Exception as e:
        st.error(f"Error reading file: {e}")
else:
    st.info("Please upload a CSV file to get started.")
    
    # Example data button
    if st.button("Load Example Data"):
        # Create a dummy dataframe for demonstration
        data = {
            "Substrate_mM": [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
            "Rate_uM_s": [0.05, 0.09, 0.22, 0.38, 0.60, 0.88, 1.05]
        }
        df_example = pd.DataFrame(data)
        csv = df_example.to_csv(index=False).encode('utf-8')
        
        st.download_button(
            label="Download Example CSV",
            data=csv,
            file_name="example_kinetics_data.csv",
            mime="text/csv",
        )
