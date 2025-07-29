// Étend l'objet global window.dash_clientside avec une nouvelle fonction personnalisée
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        // Fonction appelée lors d’un clic sur le bouton "Get SMILES"
        getKekuleSmiles: function(n_clicks) {
            // Si le bouton n'a pas encore été cliqué, on ne fait rien (Dash conserve l’état actuel)
            if (!n_clicks) return window.dash_clientside.no_update;

            // Récupère l'iframe contenant l'éditeur Kekule
            const iframe = document.getElementById("kekule-editor");

            // Vérifie que l'iframe et la fonction getSmiles sont accessibles
            if (iframe && iframe.contentWindow && typeof iframe.contentWindow.getSmiles === "function") {
                // Appelle la fonction getSmiles définie dans le fichier HTML de l’éditeur
                const smiles = iframe.contentWindow.getSmiles();
                return smiles;
            }

            // Si quelque chose ne va pas, on retourne une chaîne vide
            return '';
        }
    }
});
