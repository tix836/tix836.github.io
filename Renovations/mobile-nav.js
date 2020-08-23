let mobile = document.getElementById("_mobile-nav");
let mobileLinksButton = document.getElementById("_mobile-button-links");
let closeButton = document.getElementById("_close");
let mobileContent = document.getElementById("_mobile-nav-content");

mobileLinksButton.onclick  = function(event) {
    
    mobile.style.width = "100%";
    mobileContent.style.display = "flex"; 
}

closeButton.onclick = function(event) {
    mobile.style.width = "0%";
    mobileContent.style.display = "none";
}
