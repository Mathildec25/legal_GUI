document.addEventListener("DOMContentLoaded", function () {
  const button = document.querySelector('[id$="confirm-draw"]');
  if (button) {
    button.addEventListener("click", () => {
      const iframe = document.querySelector('iframe[id$="ketcher-frame"]');
      if (iframe && iframe.contentWindow) {
        iframe.contentWindow.postMessage("getSmiles", "*");
      }
    });
  }
});
