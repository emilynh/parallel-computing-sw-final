# Hướng Dẫn Sử Dụng SmithWaterman_GPU
## Điều Kiện Set Up Trước Khi Chạy
Để có thể chạy thành công ứng dụng, bạn cần thực hiện các bước set up sau:
Cài đặt CUDA và trình biên dịch nvcc. Xem hướng dẫn cài đặt chi tiết [Tại đây](https://drive.google.com/file/d/116fM9du8O-RDzH_G_jkNObnTDc71Z5Kt/view?usp=sharing)

## Các Bước Cài Đặt và Chạy Code
### 1. Clone Repository:
- Copy đường link repo.git: https://github.com/tklee-123/SmithWaterman_GPU.git.
- Mở terminal và chạy lệnh sau để sao chép repository về máy local:
```
git clone https://github.com/tklee-123/SmithWaterman_GPU.git
```

### 2. Mở Folder và Di Chuyển tới Thư Mục Web:
- Mở folder vừa clone về trong VSCode.
- Mở terminal trong VSCocde, di chuyển tới thư mục web:
```
cd web
```

### 3. Cài Đặt Dependencies Cho Server:
- Di chuyển vào thư mục SmithWaterman_GPU/web/server và chạy các lệnh sau:
```
cd server
npm install express
node server.js
```
- Sau khi chạy thành công lệnh inline code `node server.js`, terminal sẽ hiển thị:
> Server running on port 8000
- Tiếp theo, mở một terminal mới và di chuyển vào lại thư mục SmithWaterman_GPU/web/server
```
cd web
cd server
npm install react-scripts --save
```
> [!WARNING]
> Khi chạy inline code `npm install react-scripts --save`, nếu xảy ra lỗi tương tự như sau:
```
macbooks-MacBook-Pro-2:server macbook$ npm install react-scripts --save
npm ERR! code EEXIST
npm ERR! syscall mkdir
npm ERR! path /Users/macbook/.npm/_cacache/content-v2/sha512/fd/b7
npm ERR! errno EEXIST
npm ERR! Invalid response body while trying to fetch https://registry.npmjs.org/mini-css-extract-plugin: EACCES: permission denied, mkdir '/Users/macbook/.npm/_cacache/content-v2/sha512/fd/b7'
npm ERR! File exists: /Users/macbook/.npm/_cacache/content-v2/sha512/fd/b7
npm ERR! Remove the existing file and try again, or run npm
npm ERR! with --force to overwrite files recklessly.

npm ERR! A complete log of this run can be found in: /Users/macbook/.npm/_logs/2023-12-08T07_40_24_327Z-debug-0.log
macbooks-MacBook-Pro-2:server macbook$ 
```
- Lỗi này là do vấn đề quyền truy cập (permission denied) khi thực hiện lệnh npm install. Thông báo lỗi cụ thể là:
`EACCES: permission denied, mkdir '/Users/macbook/.npm/_cacache/content-v2/sha512/fd/b7'`
- Dưới đây là cách để giải quyết vấn đề này:

   - Chạy npm với quyền đặc quyền: Sử dụng lệnh sudo để chạy npm với quyền đặc quyền. Tuy nhiên, cần lưu ý rằng việc này có thể tạo ra các vấn đề về quyền sở hữu trong thư mục .npm. Bạn có thể chạy lệnh sau:
```
sudo npm install react-scripts --save
```

### 4. Quay Trở Lại Thư Mục Web và Cài Đặt Dependencies Cho React App:
```
cd ..
npm install
npm install react-router-dom
npm start
```

### 5. Mở Giao Diện React và Chạy Ứng Dụng:
- Khi giao diện React hiển thị, bạn sẽ thấy Home page (trang chủ) của dự án chứa các thông tin giới thiệu chung. Trên Nav bar, bạn có thể click vào "APPLICATION" hoặc trên hero page click vào nút "Get Start" để di chuyển tới phần ứng dụng.

### 6. Giao Diện Phần Ứng Dụng (APPLICATION):
- Khi giao diện phần ứng dụng hiển thị, bạn sẽ thấy thông tin giới thiệu về ứng dụng và hướng dẫn chạy.

   - Bước 1: Click vào "Choose File", chọn file test có sẵn trong folder của dự án, file test có tên inline code `a.txt` hoặc inline code `b.txt` nằm trong thư mục SmithWaterman_GPU/src hoặc tự tạo file test theo yêu cầu. File test là file chứa chuỗi ADN với các ký tự là các nucleotit.
  
   - Bước 2: Sau khi hoàn tất chọn file, click vào "Run" để chạy thuật toán. Đợi khoảng 4-5 giây để màn hình hiển thị kết quả.
  
   - Bước 3: Đối chiếu kết quả nhận được với bảng trên màn hình để xác định đặc trưng gen của chuỗi ADN truyền vào.

### 7. Xử Lý Lỗi Hiển Thị:

Nếu có lỗi hiển thị, vui lòng kiểm tra và đảm bảo rằng bạn đã cài đặt thành công môi trường để chạy thuật toán.





