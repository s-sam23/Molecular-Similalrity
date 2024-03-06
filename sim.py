from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import streamlit as st
from main import tanimoto_calc
import pandas as pd
import numpy as np



st.header('Tanimoto Moleular Similarity')
 

option = st.selectbox(
    'How would you like to be calculate the similarity',
    ('2 SMILES', 'Dataframe with SMILES column'))
if option == '2 SMILES':

    st.write('You selected:', option)
    # l1 = st.text_input('Enter the 1st SMILES')
    # l2 = st.text_input('Enter the 2nd SMILES')
    l3 = st.text_area('Enter SMILES list here : ')
    l31 = l3.split()


    if st.button('Calculate Similarity'):
        # if l1 and  l2:

        #     s = tanimoto_calc(l1,l2)
        #     st.write("Molecular Similarity for given 2 SMILES is :", s)
        if l3:
            s = tanimoto_calc(l31[0],l31[1])
            st.write("Molecular Similarity for given 2 SMILES is :", s)

if option == 'Dataframe with SMILES column':
    file = st.file_uploader('Upload your dataframe Here')
    if file:
        st.success('File uploaded successfully')
        df = pd.read_csv(file)
        st.dataframe(df)
        if st.button('Submit Dataframe'):
            # Efficient calculation using NumPy vectorization and broadcasting
            tanimoto_matrix = np.empty((len(df), len(df)))
            st.write('Shape of Dataframe  :',df.shape)
            with st.spinner():
                for i in range(len(df)):
                    mol1_smiles = df['SMILES'].iloc[i]
                    tanimoto_matrix[i] = np.array([tanimoto_calc(mol1_smiles, mol2) for mol2 in df['SMILES']])

                # Create DataFrame from the similarity matrix
                d1 = pd.DataFrame(tanimoto_matrix, columns=df['SMILES'], index=df['SMILES'])
                st.dataframe(d1)


                download_button = st.download_button(
                label="Download Data",
                data= d1.to_csv(index=False), # Convert DataFrame to CSV string
                file_name="data.csv", # Default file name
                mime="text/csv", # MIME type for CSV download
                )

                # Use an if statement to handle the button click
                if download_button:
                    st.success("Data downloaded successfully!")


        




