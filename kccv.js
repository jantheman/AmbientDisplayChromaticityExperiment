class matrix3x3 {
    constructor(m) {
        this.a = m[0][0]; this.b = m[0][1]; this.c = m[0][2];
        this.d = m[1][0]; this.e = m[1][1]; this.f = m[1][2];
        this.g = m[2][0]; this.h = m[2][1]; this.i = m[2][2];
    }

    getDeterminant() {
        return this.a * (this.e * this.i - this.f * this.h)
            - this.b * (this.d * this.i - this.f * this.g)
            + this.c * (this.d * this.h - this.e * this.g);
    }

    invert() {
        const det = this.getDeterminant();

        if (det == 0)
            throw new Error("Matrix not invertible.");

        return new matrix3x3([                          
            [(this.e * this.i - this.f * this.h) / det, (this.c * this.h - this.b * this.i) / det, (this.b * this.f - this.c * this.e) / det],
            [(this.f * this.g - this.d * this.i) / det, (this.a * this.i - this.c * this.g) / det, (this.c * this.d - this.a * this.f) / det],
            [(this.d * this.h - this.e * this.g) / det, (this.b * this.g - this.a * this.h) / det, (this.a * this.e - this.b * this.d) / det]
        ]);
    }

    transform(color) {
        let y1 = this.a * color[0] + this.b * color[1] + this.c * color[2];
        let y2 = this.d * color[0] + this.e * color[1] + this.f * color[2];
        let y3 = this.g * color[0] + this.h * color[1] + this.i * color[2];
        return [y1, y2, y3];
    }

    translateTransform(color) {
        let y1 = this.a * color[0] + this.d * color[1] + this.g * color[2];
        let y2 = this.b * color[0] + this.e * color[1] + this.h * color[2];
        let y3 = this.c * color[0] + this.f * color[1] + this.i * color[2];
        return [y1, y2, y3];
    }
}


class Kccv {
    cal_init(caldata) {
        // reshape [rY gY bY] into [rY rx ry]
        let calibration = caldata.readings.map(level => [
            [level[0][0], level[1][0], level[2][0]],
            [level[0][1], level[1][1], level[2][1]],
            [level[0][2], level[1][2], level[2][2]]
        ]);

        // average of black in XYZ
        this.baselineXYZ = this.YxyToXYZ(this.div(this.add(...calibration[0]), calibration[0].length));

        // fix data glow - for each calibration row, convert it to XYZ, subtract black of corresponding channels and convert back to Yxy
        const blackXYZ = calibration[0].map(this.YxyToXYZ);
        this.calibration = calibration = calibration.map(level => level.map((channel, index) => this.XYZtoYxy(this.sub(this.YxyToXYZ(channel), blackXYZ[index]))));
        
        // prepare expanded lookup tables
        const expandedSteps = calibration.length > 256 ? 20480 : 2048;
        this.ideal_increment = calibration[this.calibration.length - 1].map(channel => channel[0] / (expandedSteps - 1)); // channel[0] is Y
        this.gun_lut = [];
        for (var i = 0; i < expandedSteps; i++) {
            let absdiff = calibration.map(level => level.map((channel, index) => Math.abs(channel[0] - i * this.ideal_increment[index])));
            this.gun_lut[i] = this.indicesOfMin(absdiff);
        }

        // create conversion matrices
        let RGB = calibration[calibration.length - 1].map(channel => this.YxyToXYZ([1, channel[1], channel[2]]));
        this.mRGBtoXYZ = new matrix3x3(RGB);
        this.mXYZtoRGB = this.mRGBtoXYZ.invert();
    }

    add(colors) {
        let color = [...arguments[0]]; // copy
        for (let c = 1; c < arguments.length; c++)
            for (let i = 0; i < arguments[c].length; i++)
                color[i] += arguments[c][i];

        return color;
    }
    sub(colors) {
        let color = [...arguments[0]]; // copy
        for (let c = 1; c < arguments.length; c++)
            for (let i = 0; i < arguments[c].length; i++)
                color[i] -= arguments[c][i];

        return color;
    }
    mul(color, factor) {
        return color.map(c => c * factor);
    }
    div(color, factor) {
        return color.map(c => c / factor);
    }
    indicesOfMin(array) {
        let indices = new Array(array[0].length);
        for (var i = 0; i < indices.length; i++)
            indices[i] = 0;

        for (var a = 0; a < array.length; a++)
            for (var i = 0; i < indices.length; i++)
                if (array[a][i] < array[indices[i]][i])
                    indices[i] = a;

        return indices;
    }

    YxyToXYZ(Yxy) {
        const Y = Yxy[0], x = Yxy[1], y = Yxy[2];

        if (Y === 0)
            return [0, 0, 0];

        let X = x * Y / y;
        let Z = (1 - x - y) * Y / y;

        return [X, Y, Z];
    }

    XYZtoYxy(XYZ) {
        const X = XYZ[0], Y = XYZ[1], Z = XYZ[2];

        let norm = X + Y + Z;
        if (norm === 0)
            return [0, 0, 0];

        let x = X / norm;
        let y = Y / norm;

        return [Y, x, y];
    }

    LABtoXYZ(LAB, whitepointXYZ) {
        const L = LAB[0], a = LAB[1], b = LAB[2];
        const Xr = whitepointXYZ[0], Yr = whitepointXYZ[1], Zr = whitepointXYZ[2];

        const ε = 216 / 24389;
        const κ = 24389 / 27;

        let fy = (L + 16) / 116;
        let fx = a / 500 + fy;
        let fz = fy - b / 200;

        let fy3 = fy * fy * fy;
        let fx3 = fx * fx * fx;
        let fz3 = fz * fz * fz;

        let xr = fx3 > ε ? fx3 : (116 * fx - 16) / κ;
        let yr = L > κ * ε ? fy3 : L / κ;
        let zr = fz3 > ε ? fz3 : (116 * fz - 16) / κ;

        let X = xr * Xr;
        let Y = yr * Yr;
        let Z = zr * Zr;

        return [X, Y, Z];
    }

    ui8rgbToXYZ(rgb) {
        switch (this.calibration.length) {
            case 256:
                return this.uiNrgbToXYZ(rgb);
            case 1024:
                return this.uiNrgbToXYZ(this.mul(rgb, 4)); // maps 0-255 only to 0-1020, consistent with matlab kccv
            default:
                throw new Error(`Only 8-bit or 10-bit calibration data supported.`);
        }
    }
    ui10rgbToXYZ(rgb) {
        switch (this.calibration.length) {
            case 256:
                return this.uiNrgbToXYZ(this.div(rgb, 4));
            case 1024:
                return this.uiNrgbToXYZ(rgb);
            default:
                throw new Error(`Only 8-bit or 10-bit calibration data supported.`);
        }
    }

    uiNrgbToXYZ(rgb) {
        const r = rgb[0], g = rgb[1], b = rgb[2];

        let rY = this.calibration[Math.floor(r)][0].Y;
        let gY = this.calibration[Math.floor(g)][0].Y;
        let bY = this.calibration[Math.floor(b)][0].Y;

        let xyz = this.mRGBtoXYZ.translateTransform([rY, gY, bY]);
        return this.add(xyz, this.baselineXYZ);
    }

    XYZtoui8rgb(xyz, params) {
        let rgb = this.XYZtouiNrgb(xyz, params);
        switch (this.calibration.length) {
            case 256:
                return rgb;
            case 1024:
                return rgb.map(v => Math.round(v / 4));
            default:
                throw new Error(`Only 8-bit or 10-bit calibration data supported.`);
          }
    }

    XYZtoui10rgb(xyz, params) {
        let rgb = this.XYZtouiNrgb(xyz, params);
        switch (this.calibration.length) {
            case 256:
                return this.mul(rgb, 4);
            case 1024:
                return rgb;
            default:
                throw new Error(`Only 8-bit or 10-bit calibration data supported.`);
        }
    }

    XYZtouiNrgb(xyz, params) {
        xyz = this.sub(xyz, this.baselineXYZ);
        let rgb = this.kccv(xyz, 'XYZ', 'RGB');

        if (params?.spat_correct) {
            throw new Error('Correction of spatial non-uniformaities not implemented.');
        }

        rgb = rgb.map((v, index) => v / this.ideal_increment[index]);

        if (params?.RGB_Ceil) {
            // remove the undisplayable and replace with either 0 or ceiling values
            rgb = rgb.map(v => Math.min(2047, Math.max(0, v)));
        } else {
            // out of gamut
            if (rgb.some(v => v < 0 || v > 2047)) return [0, 0, 0];
        }

        // round the entries so they can serve as indexes to a matrix
        rgb = rgb.map(v => Math.floor(v));

        // remove NaN values (little x or y values == 0)
        // if (rgb.some(v => Number.isNaN(v)) return [0, 0, 0];

        if (params?.dither) {
            throw new Error('Dithering not implemented.');
        }
        else {
            rgb = rgb.map((v, index) => this.gun_lut[v][index]);
        }

        return rgb;
    }

    toXYZ(color, from, params) {
        switch (from.toLowerCase()) {
            case 'rgb': return this.mRGBtoXYZ.translateTransform(from);
            case 'yxy': return this.YxyToXYZ(color);
            case 'ui8rgb': return this.ui8rgbToXYZ(color);
            case 'ui10rgb': return this.ui10rgbToXYZ(color);
            case 'lab': return this.LABtoXYZ(color, this.YxyToXYZ(params.whitepoint));
            case 'xyz': return [...color]; // copy
            default:
                throw new Error(`Conversion from '${from}' not implemented.`);
        }
    }

    kccv(color, from, to, params) {
        let xyz = this.toXYZ(color, from, params);
        switch (to.toLowerCase()) {
            case 'xyz': return xyz;
            case 'rgb': return this.mXYZtoRGB.translateTransform(xyz);
            case 'yxy': return this.XYZtoYxy(xyz);
            case 'ui8rgb': return this.XYZtoui8rgb(xyz, params);
            case 'ui10rgb': return this.XYZtoui10rgb(xyz, params);
            default:
                throw new Error(`Conversion to '${to}' not implemented.`);
        }
    }
}