// Extend the global window.dash_clientside object with a custom function
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        // Function triggered when the "Get SMILES" button is clicked
        getKekuleSmiles: function(n_clicks) {
            // If the button hasn't been clicked yet, do nothing (Dash keeps current state)
            if (!n_clicks) return window.dash_clientside.no_update;

            // Get the iframe containing the Kekule editor
            const iframe = document.getElementById("kekule-editor");

            // Check if the iframe and the getSmiles function are accessible
            if (iframe && iframe.contentWindow && typeof iframe.contentWindow.getSmiles === "function") {
                // Call the getSmiles function defined in the editor’s HTML file
                const smiles = iframe.contentWindow.getSmiles();
                return smiles;
            }

            // If something goes wrong, return an empty string
            return '';
        }
    }
});
