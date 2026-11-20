import Foundation
import CoreImage

let args = CommandLine.arguments
let msg = args[1]
let level = args.count > 2 ? args[2] : "Q"
let f = CIFilter(name: "CIQRCodeGenerator")!
f.setValue(msg.data(using: .utf8)!, forKey: "inputMessage")
f.setValue(level, forKey: "inputCorrectionLevel")
let img = f.outputImage!
let e = img.extent
let w = Int(e.width), h = Int(e.height)
let ctx = CIContext(options: nil)
var buf = [UInt8](repeating: 0, count: w * h * 4)
ctx.render(img, toBitmap: &buf, rowBytes: w * 4, bounds: e, format: .RGBA8, colorSpace: CGColorSpaceCreateDeviceRGB())
print("size \(w) \(h)")
for y in 0..<h {
    var line = ""
    for x in 0..<w { line += buf[(y * w + x) * 4] < 128 ? "1" : "0" }
    print(line)
}
