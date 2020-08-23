// или сделать через addEventListener


let watch = document.getElementById('_watch');
let video = document.getElementById('video_con');
let header = document.getElementById('_header');

let coords = header.getBoundingClientRect();


watch.onclick = function(event) {

    video.style.display = `block`;

      window.scrollTo({
      top: coords.bottom - 50,
      behavior: "smooth"
     });
}